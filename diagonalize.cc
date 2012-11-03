/*
 * diagonalize.cc
 */

#include "llvm/Pass.h"
#include "llvm/Function.h"
#include "llvm/Constants.h"
#include "llvm/Instructions.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/ADT/ValueMap.h"
#include "llvm/Support/PatternMatch.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <iostream>
#include "llvm/Support/raw_ostream.h"

using namespace llvm;
using namespace PatternMatch;
using namespace Eigen;

#define OP_IN_RANGE(_op, _start, _end) \
    (_op >= Instruction:: _start && _op <= Instruction:: _end)

namespace
{

struct ADPass : public LoopPass
{
private:
    Loop* loop;
    DominatorTree* DT;
    BasicBlock* exit_block;
    std::vector<BasicBlock*> blocks;

    /* Keep track of the loop range. */
    Value* iter_final;
    PHINode* iter_var;

    /* Index and store state variables. */
    ValueMap<PHINode*, size_t> phis;

    /* Track linear combinations of variables. */
    typedef ValueMap<PHINode*, double> Coefficients;

public:
    static char ID;

    ADPass() : LoopPass(ID) {}

    virtual void getAnalysisUsage(AnalysisUsage& AU) const {
        AU.addRequired<DominatorTree>();
    }

    virtual bool runOnLoop(Loop* L, LPPassManager&) {
        loop = L;
        blocks = loop->getBlocks();
        DT = &getAnalysis<DominatorTree>();
        phis.clear();

        if (!loop->isLCSSAForm(*DT)
            || loop->getSubLoops().size()
            || !(exit_block = loop->getUniqueExitBlock()))
        {
            return false;
        }

        /* Filter away loops with unsupported instructions. */
        size_t nr_cmps = 0;
        ICmpInst* icmp = NULL;
        for (auto it = blocks.begin(); it != blocks.end(); ++it) {
            BasicBlock* BB = *it;
            for (auto II = BB->begin(); II != BB->end(); ++II) {
                Instruction* instr = II;
                if (isa<ICmpInst>(instr)) {
                    if (++nr_cmps > 1) {
                        return false;
                    } else {
                        icmp = cast<ICmpInst>(instr);
                    }
                } else if (isa<PHINode>(instr)) {
                    PHINode* PN = cast<PHINode>(instr);
                    if (PN->getNumIncomingValues() != 2
                        || !DT->properlyDominates(PN->getIncomingBlock(0),
                                                  blocks.front()))
                    {
                        return false;
                    }
                    Value* inLhs = PN->getIncomingValue(0);
                    if (!(isa<ConstantInt>(inLhs) || isa<ConstantFP>(inLhs))) {
                        return false;
                    } else {
                        phis[PN] = 0;
                    }
                } else if (isa<BinaryOperator>(instr)) {
                    BinaryOperator* binop = cast<BinaryOperator>(instr);
                    if (!(OP_IN_RANGE(binop->getOpcode(), Add, FDiv))) {
                        return false;
                    }
                } else if (!(isa<BranchInst>(instr)
                            || isa<IndirectBrInst>(instr)))
                {
                    return false;
                }
            }
        }

        /* Extract the loop iterator if possible. */
        if (!icmp || !extractIterator(icmp)) {
            return false;
        } else {
            phis.erase(iter_var);
        }

        size_t phi_index = 0;
        for (auto kv = phis.begin(); kv != phis.end(); ++kv) {
            phis[kv->first] = phi_index++;
        }

        /* Find the initial state and all linear dependence relations. */
        size_t nr_phis = phis.size();
        MatrixXd InitialState(nr_phis, 1);
        MatrixXd TransformationMatrix(nr_phis, nr_phis);
        TransformationMatrix << MatrixXd::Zero(nr_phis, nr_phis);
        for (auto kv = phis.begin(); kv != phis.end(); ++kv) {
            PHINode* PN = kv->first;
            size_t phi_label = phis[PN];
            InitialState(phi_label, 0) = ToDouble(PN->getIncomingValue(0));

            Coefficients coeffs;
            if (!trackUpdates(PN->getIncomingValue(1), coeffs)) {
                return false;
            }
            for (auto ckv = coeffs.begin(); ckv != coeffs.end(); ++ckv) {
                size_t target_phi = phis[ckv->first];
                TransformationMatrix(phi_label, target_phi) = ckv->second;
            }
        }

        /* Diagonalize the transformation matrix. */
        EigenSolver<MatrixXd> EigSolver(TransformationMatrix);
        MatrixXcd P = EigSolver.eigenvectors();
        MatrixXcd D = EigSolver.eigenvalues().asDiagonal();
        MatrixXcd Pinv = P.inverse();
        if (!checkSystem(TransformationMatrix, InitialState, P, D, Pinv)) {
            return false;
        }

        /* Emit instructions to compute the closed form in a new block. */
        LLVMContext& ctx = blocks.front()->getContext();
        Module* mod = exit_block->getParent()->getParent();
        BasicBlock* dgen = BasicBlock::Create(ctx, "dgen", 0, exit_block);
        Type* numTy = Type::getDoubleTy(ctx);
        std::vector<Type*> powProto(2, numTy);
        FunctionType* powType = FunctionType::get(numTy, powProto, false);
        Function* powf = Function::Create(powType,
            GlobalValue::ExternalLinkage, "llvm.pow.f64", mod);
        powf->setCallingConv(CallingConv::C);

        Value** PDn = new Value*[nr_phis * nr_phis];
        for (size_t j = 0; j < nr_phis; ++j) {
            Value* exptargs[] = {
                ToConstantFP(ctx, std::real(D(j, j))), iter_final
            };
            Value* eigvexpt = CallInst::Create(powf,
                ArrayRef<Value*>(exptargs, 2), "eigvexpt", dgen);

            for (size_t i = 0; i < nr_phis; ++i) {
                size_t index = i * nr_phis + j;
                PDn[index] = BinaryOperator::Create(Instruction::FMul,
                    ToConstantFP(ctx, std::real(P(i, j))), eigvexpt, "pdn",
                    dgen);
            }
        }

        Value** soln = new Value*[nr_phis];
        for (size_t i = 0; i < nr_phis; ++i) {
            for (size_t j = 0; j < nr_phis; ++j) {
                Value* inprod = ToConstantFP(ctx, 0.0);
                for (size_t k = 0; k < nr_phis; ++k) {
                    Value* ik_kj = BinaryOperator::Create(Instruction::FMul,
                        PDn[i * nr_phis + k],
                        ToConstantFP(ctx, std::real(Pinv(k, j))),
                        "ik_kj", dgen);
                    inprod = BinaryOperator::Create(Instruction::FAdd,
                        ik_kj, inprod, "inprod", dgen);
                }

                if (soln[i] == NULL) {
                    soln[i] = inprod;
                } else {
                    Value* xj_prod = BinaryOperator::Create(Instruction::FMul,
                        inprod, ToConstantFP(ctx, InitialState(j)), "pdpxj",
                        dgen);
                    soln[i] = BinaryOperator::Create(Instruction::FAdd, 
                        xj_prod, soln[i], "xf", dgen);
                }
            }
        }

        errs() << dgen << "\n";

        delete[] PDn;
        delete[] soln;

        /*
         * XXX:
         * - Do some analysis on the loop;
         *   1) What is the block we jump out to when we finish the loop?
         *      exit_block
         *   2) Which PHI nodes in this block have incoming values from
         *      our loop? Store these PHI nodes: determine the state variable
         *      they refer to by examining the Value* they have incoming, and
         *      searching for that Value* in our set of PHIs.
         *      Map[OuterPHI] => StateVariable
         *   3) Delete all of the Instructions in our loop.
         *   4) Inject matrix multiplication code into the hollowed BB.
         *   5) Change the incoming values in the OuterPHIs to the final
         *      results produced by the diagonalization process.
         */

        /* Find dependencies on state variables outside of the loop. */
        ValueMap<PHINode*, PHINode*> outerDeps;
        for (auto II = exit_block->begin(); II != exit_block->end(); ++II) {
            Instruction* instr = II;
            if (isa<PHINode>(instr)) {
                PHINode* PN = cast<PHINode>(instr);
                if (PN->getBasicBlockIndex(blocks.back()) != -1) {
                    Value* V = PN->getIncomingValueForBlock(blocks.back());
                    for (auto kv = phis.begin(); kv != phis.end(); ++kv) {
                        if (kv->first->getIncomingValue(1) == V) {
                            outerDeps[PN] = kv->first;
                        }
                    }
                }
            }
        }

        return false;
    }

private:
    bool extractIterator(ICmpInst* icmp) {
        /* In canonical form, we should get an icmp with the predicate
         * reduced to an equality test. */
        Value *cmpLhs, *cmpRhs;
        ICmpInst::Predicate IPred;
        if (match(icmp, m_ICmp(IPred, m_Value(cmpLhs),
                                      m_Value(cmpRhs)))
            && IPred == CmpInst::Predicate::ICMP_EQ)
        {
            /* We should be comparing against the final iterator value, which
             * should not be found in the loop body. */
            if (isa<Instruction>(cmpRhs)) {
                Instruction* instr = cast<Instruction>(cmpRhs);
                BasicBlock* BB = instr->getParent();
                if (!DT->properlyDominates(BB, blocks.front())) {
                    return false;
                }
                if (isa<Constant>(instr->getOperand(0))
                    || isa<Instruction>(instr->getOperand(0)))
                {
                    return false;
                }
                iter_final = instr->getOperand(0);
            } else {
                iter_final = cmpRhs;
            }

            /* Check if the loop increment is behaving as we expect it to. */
            Value *incrLhs, *incrRhs;
            if (match(cmpLhs, m_Add(m_Value(incrLhs), m_Value(incrRhs)))) {
                if (isa<PHINode>(incrLhs)
                    && isa<ConstantInt>(incrRhs)
                    && cast<ConstantInt>(incrRhs)->isOne())
                {
                    iter_var = cast<PHINode>(incrLhs);
                    Value* iter_initial = iter_var->getIncomingValue(0);
                    if (isa<ConstantInt>(iter_initial)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    bool trackUpdates(Value* parent, Coefficients& coeffs) {
        /* Determine the linear combination that produces 'parent'. */
        if (isa<Constant>(parent)) {
            return true;
        } else if (isa<PHINode>(parent)) {
            /* A PHINode should contribute a single copy of itself. */
            PHINode* PN = cast<PHINode>(parent);
            if (!phis.count(PN)) {
                return false;
            } else {
                coeffs[PN] = 1;
            }
        } else if (isa<BinaryOperator>(parent)) {
            BinaryOperator* binop = cast<BinaryOperator>(parent);
            int opcode = binop->getOpcode();
            Value *LHS = binop->getOperand(0), *RHS = binop->getOperand(1);

            Coefficients lhsCoeffs, rhsCoeffs;
            if (!trackUpdates(LHS, lhsCoeffs)
                || !trackUpdates(RHS, rhsCoeffs))
            {
                return false;
            }

            double scalar;
            if (OP_IN_RANGE(opcode, Add, FSub)) {
                /* Add instructions shouldn't operate on scalars. */
                if (isScalar(lhsCoeffs) || isScalar(rhsCoeffs)) {
                    return false;
                }
            } else {
                /* Mul instructions can only have one scalar operand. */
                if (isScalar(lhsCoeffs) ^ isScalar(rhsCoeffs)) {
                    return false;
                }

                /* Div instructions cannot have scalar numerators. */
                if (OP_IN_RANGE(opcode, UDiv, FDiv) && !isScalar(lhsCoeffs)) {
                    return false;
                }

                scalar = ToDouble(isScalar(lhsCoeffs) ? LHS : RHS);
            }

            /* Merge the two sets of coefficients. */
            for (auto kv = phis.begin(); kv != phis.end(); ++kv) {
                PHINode* PN = kv->first;
                double lcoeff = lhsCoeffs.lookup(PN);
                double rcoeff = rhsCoeffs.lookup(PN);

                /* Adding nil entries to 'coeffs' breaks isScalar(X). */
                if (lcoeff == 0.0 && rcoeff == 0.0) {
                    continue;
                }

                if (OP_IN_RANGE(opcode, Add, FAdd)) {
                    coeffs[PN] = lcoeff + rcoeff;
                } else if (OP_IN_RANGE(opcode, Sub, FSub)) {
                    coeffs[PN] = lcoeff - rcoeff;
                } else if (OP_IN_RANGE(opcode, Mul, FMul)) {
                    coeffs[PN] = scalar * (lcoeff + rcoeff);
                } else {
                    coeffs[PN] = (1 / scalar) * (lcoeff + rcoeff);
                }
            }
        } else return false;
        return true;
    }

    int64_t ToInt(Value* V) {
        return cast<ConstantInt>(V)->getValue().getSExtValue();
    }

    double ToDouble(Value* V) {
        if (isa<ConstantInt>(V)) {
            return double(ToInt(V));
        } else if (isa<ConstantFP>(V)) {
            return cast<ConstantFP>(V)->getValueAPF().convertToDouble();
        }
        return 0.0;
    }

    Value* ToConstantFP(LLVMContext& ctx, double n) {
        return ConstantFP::get(ctx, APFloat(n));
    }

    bool isScalar(Coefficients& coeffs) {
        return coeffs.size() == 0;
    }

    bool checkSystem(MatrixXd& A, MatrixXd& x,
                     MatrixXcd& P, MatrixXcd& D, MatrixXcd& Pinv)
    {
        /* Ensure that there are no complex eigenvectors. */
        for (int i = 0; i < P.rows(); ++i) {
            for (int j = 0; j < P.cols(); ++j) {
                if (std::imag(P(i, j)) != 0.0) {
                    return false;
                }
            }
        }

        /* Check that the diagonalization is reasonably accurate. */
        MatrixXd A3x = A*A*A*x;
        MatrixXcd PD3Px = (P*(D*D*D)*Pinv)*x;
        const double epsilon = 100 * std::numeric_limits<double>::epsilon();
        for (int i = 0; i < x.rows(); ++i) {
            if (epsilon <
                std::abs(std::abs(PD3Px(i, 0)) - std::abs(A3x(i, 0))))
            {
                return false;
            }
        }
        return true;
    }
};

char ADPass::ID = 0;

}

static RegisterPass<ADPass> X("auto-diagonalize",
    "Diagonalize linear dynamical systems", false, false);

/*
 * diagonalize.cc
 */

#include "llvm/Pass.h"
#include "llvm/Function.h"
#include "llvm/Module.h"
#include "llvm/Constants.h"
#include "llvm/Instructions.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/ADT/ValueMap.h"
#include "llvm/Support/PatternMatch.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

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
    ICmpInst* icmp; 
    Value* iter_final;
    ConstantInt* iter_off;
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
            || !(exit_block = loop->getUniqueExitBlock())
            || !loopFilter())
        {
            return false;
        }

        /* Extract the loop iterator if possible. */
        if (!icmp || !extractIterator()) {
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
        if (!checkSystem(TransformationMatrix, P, D, Pinv)) {
            return false;
        }

        /* Emit instructions to compute the closed form in a new block. */
        LLVMContext& ctx = blocks.front()->getContext();
        Function* parentFunc = exit_block->getParent();
        Module* mod = parentFunc->getParent();
        BasicBlock* dgen = BasicBlock::Create(ctx, "dgen", parentFunc,
            exit_block);
        Type* numTy = Type::getDoubleTy(ctx);
        Function* powf = NULL;
        if (!(powf = mod->getFunction("llvm.pow.f64"))) {
            std::vector<Type*> powProto(2, numTy);
            FunctionType* powType = FunctionType::get(numTy, powProto, false);
            powf = Function::Create(powType, GlobalValue::ExternalLinkage,
                "llvm.pow.f64", mod);
            powf->setCallingConv(CallingConv::C);
        }

        /* P(D^n) = r * λ^n : ∀(r) ∈ row(P) */
        Value** PDn = new Value*[nr_phis * nr_phis];
        Value* iexpt = BinaryOperator::Create(Instruction::Sub, iter_final,
           iter_off, "iexpt", dgen);
        Value* exponent = CastInst::Create(Instruction::UIToFP, iexpt, numTy,
            "fexpt", dgen);
        for (size_t j = 0; j < nr_phis; ++j) {
            /* λ[j]^n */
            Value* exptargs[] = {
                ToConstantFP(ctx, std::real(D(j, j))), exponent
            };
            Value* eigvexpt = CallInst::Create(powf,
                ArrayRef<Value*>(exptargs, 2), "eigvexpt", dgen);

            /* PDn[i][j] = P[i][j] * λ[j]^n */
            for (size_t i = 0; i < nr_phis; ++i) {
                size_t index = i * nr_phis + j;
                PDn[index] = BinaryOperator::Create(Instruction::FMul,
                    ToConstantFP(ctx, std::real(P(i, j))), eigvexpt, "pdn",
                    dgen);
            }
        }

        /* xf = P(D^n) * Pinv * x0 */
        Value** soln = new Value*[nr_phis];
        Value* zero = ToConstantFP(ctx, 0.0);
        for (size_t i = 0; i < nr_phis; ++i) {
            soln[i] = zero;
            for (size_t j = 0; j < nr_phis; ++j) {
                /* dotp = <a, b> : a ∈ row(P(D^n)), b ∈ col(Pinv) */
                Value* dotp = zero;
                for (size_t k = 0; k < nr_phis; ++k) {
                    Value* ik_kj = BinaryOperator::Create(Instruction::FMul,
                        PDn[i * nr_phis + k],
                        ToConstantFP(ctx, std::real(Pinv(k, j))),
                        "ik_kj", dgen);
                    dotp = BinaryOperator::Create(Instruction::FAdd,
                        ik_kj, dotp, "dotp", dgen);
                }

                /* xf[i] = ∑ P(D^n)Pinv[i][j] * x0[j] */
                Value* xj_prod = BinaryOperator::Create(Instruction::FMul,
                    dotp, ToConstantFP(ctx, InitialState(j)), "pdpxj",
                    dgen);
                soln[i] = BinaryOperator::Create(Instruction::FAdd, 
                    xj_prod, soln[i], "xf", dgen);
            }
        }

        delete[] PDn;

        /* Point outgoing edges from the loop preheader to dgen. */
        BranchInst::Create(exit_block, dgen);
        BasicBlock* preheader = loop->getLoopPreheader();
        TerminatorInst* TI = preheader->getTerminator();
        BranchInst* br = dyn_cast<BranchInst>(TI);
        br->setSuccessor(0, dgen);

        /* Replace dependencies on state variables outside of the loop. */
        for (auto II = exit_block->begin(); II != exit_block->end(); ++II) {
            Instruction* instr = II;
            if (isa<PHINode>(instr)) {
                PHINode* exitPhi = cast<PHINode>(instr);
                int incomingEdge = exitPhi->getBasicBlockIndex(blocks.back());
                if (incomingEdge != -1) {
                    Value* exitV = exitPhi->getIncomingValue(incomingEdge);
                    for (auto kv = phis.begin(); kv != phis.end(); ++kv) {
                        PHINode* loopPhi = kv->first;
                        if (loopPhi->getIncomingValue(1) != exitV) {
                            continue;
                        }
                        size_t phiIdx = phis[loopPhi];
                        exitPhi->setIncomingValue(incomingEdge, soln[phiIdx]);
                        exitPhi->setIncomingBlock(incomingEdge, dgen);
                    }
                }
            }
        }

        /* Delete the old loop. */
        for (auto it = blocks.begin(); it != blocks.end(); ++it) {
            BasicBlock* BB = *it;
            BB->eraseFromParent();
        }

        delete[] soln;

        return true;
    }

private:
    bool loopFilter() {
        /* Filter away loops with unsupported instructions. */
        icmp = NULL;
        size_t nr_cmps = 0;
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
                } else if (!isa<BranchInst>(instr)) {
                    return false;
                }
            }
        }
        return true;
    }

    bool extractIterator() {
        /* In canonical form, we should get an icmp with the predicate
         * reduced to an equality test. */
        Value* loop_var;
        ConstantInt* loop_incr;
        ICmpInst::Predicate IPred;
        if (!match(icmp, m_ICmp(IPred, m_Add(m_Value(loop_var),
                                             m_ConstantInt(loop_incr)),
                                       m_Value(iter_final)))
            || IPred != CmpInst::Predicate::ICMP_EQ
            || !isa<PHINode>(loop_var)
            || !loop_incr->isOne())
        {
            return false;
        }

        /* We need to find the loop offset to determine our final exponent. */
        iter_var = cast<PHINode>(loop_var);
        if (!match(iter_var->getIncomingValue(0), m_ConstantInt(iter_off))) {
            return false;
        }
        return true;
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
                if (!(isScalar(lhsCoeffs) ^ isScalar(rhsCoeffs))) {
                    return false;
                }

                /* Div instructions cannot have scalar numerators. */
                if (OP_IN_RANGE(opcode, UDiv, FDiv) && isScalar(lhsCoeffs)) {
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
        } else {
            return false;
        }
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

    bool checkSystem(MatrixXd& A, MatrixXcd& P,
                     MatrixXcd& D, MatrixXcd& Pinv)
    {
        /* Check if the diagonalization worked. */
        MatrixXcd PDPi = P*D*Pinv;
        const double epsilon = 25 * std::numeric_limits<double>::epsilon();
        for (int i = 0; i < P.rows(); ++i) {
            for (int j = 0; j < P.cols(); ++j) {
                if (std::imag(P(i, j)) != 0.0) {
                    return false;
                }
                if (std::abs(A(i, j) - PDPi(i, j)) . epsilon) {
                    return false;
                }
            }
        }
        return true;
    }
};

char ADPass::ID = 0;

}

static RegisterPass<ADPass> X("auto-diagonalize",
    "Diagonalize linear dynamical systems", false, false);

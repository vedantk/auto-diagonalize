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
    std::vector<BasicBlock*> blocks;

    /* Keep track of the loop range. */
    int64_t iter_base;
    Value* iter_final;
    PHINode* iter_var;

    /* Index and store state variables. */
    ValueMap<PHINode*, int> phis;

    /* Track linear combinations of variables. */
    typedef ValueMap<PHINode*, double> Coefficients;

public:
    static char ID;

    ADPass() : LoopPass(ID) {}

    virtual void getAnalysisUsage(AnalysisUsage& AU) const {
        AU.addRequired<DominatorTree>();
    }

    virtual bool runOnLoop(Loop* L, LPPassManager&) {
        if (L->getSubLoops().size()) {
            return false;
        }

        if (!L->getUniqueExitBlock()) {
            return false;
        }

        loop = L;
        blocks = loop->getBlocks();
        DT = &getAnalysis<DominatorTree>();
        phis.clear();
        
        /* Filter away loops with unsupported instructions. */
        int nr_cmps = 0;
        ICmpInst* icmp = NULL;
        for (auto it = blocks.begin(); it != blocks.end(); ++it) {
            BasicBlock* BB = *it;
            for (auto II = BB->begin(); II != BB->end(); ++II) {
                Instruction* instr = static_cast<Instruction*>(II);
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

        int phi_index = 0;
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
            int phi_label = phis[PN];
            InitialState(phi_label, 0) = ToDouble(PN->getIncomingValue(0));

            Coefficients coeffs;
            if (!trackUpdates(PN->getIncomingValue(1), coeffs)) {
                return false;
            }
            for (auto ckv = coeffs.begin(); ckv != coeffs.end(); ++ckv) {
                int target_phi = phis[ckv->first];
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

        std::cout << "System valid;\n";
        std::cout << TransformationMatrix << std::endl;

        return false;
    }

private:
    bool extractIterator(ICmpInst* icmp) {
        /* In canonical form, LLVM should pass us an ICmp with the predicate
         * reduced to an equality comparison. The LHS should contain the loop
         * increment, and the RHS should be an out-of-loop value. */
        Value *cmpLhs, *cmpRhs;
        ICmpInst::Predicate IPred; 
        if (match(icmp, m_ICmp(IPred, m_Value(cmpLhs),
                                      m_Value(cmpRhs)))
            && IPred == CmpInst::Predicate::ICMP_EQ)
        {
            if (isa<Instruction>(cmpRhs)) {
                BasicBlock* BB = cast<Instruction>(cmpRhs)->getParent();
                if (!DT->properlyDominates(BB, blocks.front())) {
                    return false;
                }
            }

            Value *incrLhs, *incrRhs;
            if (match(cmpLhs, m_Add(m_Value(incrLhs), m_Value(incrRhs)))) {
                if (isa<PHINode>(incrLhs)
                    && isa<ConstantInt>(incrRhs)
                    && cast<ConstantInt>(incrRhs)->isOne())
                {
                    iter_final = cmpRhs;
                    iter_var = cast<PHINode>(incrLhs);
                    Value* iter_initial = iter_var->getIncomingValue(0);
                    if (isa<ConstantInt>(iter_initial)) {
                        iter_base = ToInt(iter_initial);
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
                if (scalarp(lhsCoeffs) || scalarp(rhsCoeffs)) {
                    return false;
                }
            } else {
                /* Mul instructions can only have one scalar operand. */
                if (scalarp(lhsCoeffs) ^ scalarp(rhsCoeffs)) {
                    return false;
                }

                /* Div instructions can not have scalar numerators. */
                if (OP_IN_RANGE(opcode, UDiv, FDiv) && !scalarp(lhsCoeffs)) {
                    return false;
                }

                scalar = ToDouble(scalarp(lhsCoeffs) ? LHS : RHS);
            }

            /* Merge the two sets of coefficients. */
            for (auto kv = phis.begin(); kv != phis.end(); ++kv) {
                PHINode* PN = kv->first;
                double lcoeff = lhsCoeffs.lookup(PN);
                double rcoeff = rhsCoeffs.lookup(PN);
                
                /* Adding nil entries to 'coeffs' breaks scalarp(X). */
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

    bool scalarp(Coefficients& coeffs) {
        return coeffs.size() == 0;
    }

    bool checkSystem(MatrixXd& A, MatrixXd& x,
                     MatrixXcd& P, MatrixXcd& D, MatrixXcd& Pinv)
    {
        /* Ensure that there are no complex eigenvectors. */
        for (int i=0; i < P.rows(); ++i) {
            for (int j=0; j < P.cols(); ++j) {
                if (std::imag(P(i, j)) != 0.0) {
                    return false;
                }
            }
        }

        /* Check that the diagonalization is reasonably accurate. */
        MatrixXd A3x = A*A*A*x;
        MatrixXcd PD3Px = (P*(D*D*D)*Pinv)*x;
        const double epsilon = 100 * std::numeric_limits<double>::epsilon();
        for (int i=0; i < x.rows(); ++i) {
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

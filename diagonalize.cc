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
#include "llvm/Support/raw_ostream.h"

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace llvm;
using namespace PatternMatch;
using namespace Eigen;

#define OP_IN_RANGE(_op, _start, _end) \
    (_op >= Instruction::_start && _op <= Instruction::_end)

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

    /* Store tags for state variables. */
    ValueMap<PHINode*, int> phis;
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
            errs() << "PHINode: " << *kv->first
                   << "\n\t Label: " << phi_index << "\n";
            phis[kv->first] = phi_index++;
        }
        
        /* Find the initial state and all linear dependence relations. */
        MatrixXd InitialState(phis.size(), 1);
        MatrixXcd TransformationMatrix(phis.size(), phis.size());
        for (auto kv = phis.begin(); kv != phis.end(); ++kv) {
            PHINode* PN = kv->first;
            int phi_label = phis[PN];
            InitialState(phi_label, 0) = DoubleFromValue(PN->getIncomingValue(0));

            Coefficients coeffs;
            if (!trackUpdates(phi_label, PN->getIncomingValue(1), coeffs)) {
                return false;
            }
            for (auto ckv = coeffs.begin(); ckv != coeffs.end(); ++ckv) {
                int target_phi = phis[kv->first];
                TransformationMatrix(phi_label, target_phi) = kv->second;
            }
        }

#if 0
        /* XXX:
         * How do we check if complex eigenvalues were generated? Avoid this case...
         */
        EigenSolver<MatrixXd> EigGen(TransformationMatrix);
        MatrixXcd P = EigGen.eigenvectors();
        MatrixXcd D = EigGen.eigenvalues().asDiagonal();
        MatrixXcd Pinv = P.inverse();

        /* XXX:
         * Perform P * (D^[iter_final]) * Pinv * InitialState, and assign the final
         * values back into the state variables.
         */
#endif
        
        return false;
    }

private:
    bool extractIterator(ICmpInst* icmp) {
        /* In canonical form, LLVM should pass us an ICmp with the predicate
         * reduced to an equality comparison. The LHS should contain the loop
         * increment, and the RHS should be an out-of-loop 'constant'. */
        Value *cmpLhs, *cmpRhs;
        ICmpInst::Predicate IPred; 
        if (match(icmp, m_ICmp(IPred, m_Value(cmpLhs),
                                      m_Value(cmpRhs)))
            && IPred == CmpInst::Predicate::ICMP_EQ)
        {
            if (isa<Instruction>(cmpRhs)) {
                Instruction* irhs = cast<Instruction>(cmpRhs);
                if (!DT->properlyDominates(irhs->getParent(), blocks.front())) {
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
                        iter_base = IntFromValue(iter_initial);
                        return true;
                    }
                }
            }
        }
        return false;
    }

    bool trackUpdates(Value* parent, Coefficients& coeffs) {
        if (isa<PHINode>(parent)) {
            /* Isolated PHINodes contribute one copy of themselves. */
            PHINode* PN = cast<PHINode>(parent);
            coeffs[PN] = 1;
        } else if (isa<BinaryOperator>(parent)) {
            BinaryOperator* binop = cast<BinaryOperator>(parent);
            int opcode = binop->getOpcode();
            Value *LHS = binop->getOperand(0), *RHS = binop->getOperand(1);
            Coefficients lhsCoeffs, rhsCoeffs;
            double scalar;

            if (!(trackUpdates(LHS, lhsCoeffs) && trackUpdates(RHS, rhsCoeffs))) {
                return false;
            }

            if (OP_IN_RANGE(opcode, Add, Fsub)) {
                /* Add instructions shouldn't operate on scalars. */
                if (scalarExpr(lhsCoeffs) || scalarExpr(rhsCoeffs)) {
                    return false;
                }
            } else {
                /* Multiply instructions must have (only) one scalar operand. */
                if (scalarExpr(lhsCoeffs) ^ scalarExpr(rhsCoeffs)) {
                    return false;
                }

                /* Divide instructions cannot have scalar numerators. */
                if (OP_IN_RANGE(opcode, Udiv, FDiv) && !scalarExpr(LHS)) {
                    return false;
                }

                scalar = DoubleFromValue(scalarExpr(lhsCoeffs) ? LHS : RHS);
            }

            for (auto kv = phis.begin(); kv != phis.end(); ++kv) {
                PHINode* PN = kv->first;
                double lcoeff = lhsCoeffs.lookup(PN), rcoeff = rhsCoeffs.lookup(PN);
                if (lcoeff == 0 && rcoeff == 0) {
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

    int64_t IntFromValue(Value* V) {
        return cast<ConstantInt>(V)->getValue().getSExtValue();
    }

    double DoubleFromValue(Value* V) {
        if (isa<ConstantInt>(V)) {
            return double(IntFromValue(V));
        } else if (isa<ConstantFP>(V)) {
            return cast<ConstantFP>(V)->getValueAPF().convertToDouble();
        }
    }

    bool scalarExpr(Coefficients& coeffs) {
        return coeffs.size() == 0;
    }
};

char ADPass::ID = 0;

}
 
static RegisterPass<ADPass> X("auto-diagonalize",
    "Diagonalize linear dynamical systems", false, false);

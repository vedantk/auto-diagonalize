/*
 * diagonalize.cc
 */

#include "llvm/Pass.h"
#include "llvm/Function.h"
#include "llvm/Constants.h"
#include "llvm/Instructions.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/Support/PatternMatch.h"
#include "llvm/Support/raw_ostream.h"

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace llvm;
using namespace PatternMatch;
using namespace Eigen;

namespace
{

struct ADPass : public LoopPass
{
    Loop* loop;
    DominatorTree* DT;
    std::vector<BasicBlock*> blocks;

    /* Track the loop range, iterator, and transformation matrix. */
    int64_t iter_base;
    Value* iter_final;
    PHINode* iter_var;
    MatrixXcd TransformationMatrix;

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
        
        /* Filter away loops with unsupported instructions. */
        int nr_cmps = 0;
        ICmpInst* icmp = NULL;
        SmallPtrSet<PHINode*, 4> phis;
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
                    if (!isa<Constant>(inLhs)) {
                        return false;
                    } else {
                        phis.insert(PN);
                    }
                } else if (isa<BinaryOperator>(instr)) {
                    BinaryOperator* binop = cast<BinaryOperator>(instr);
                    if (!(binop->getOpcode() >= Instruction::Add
                        && binop->getOpcode() <= Instruction::FDiv))
                    {
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
        
        /* Find the initial state and all linear dependence relations. */
        int count = 0;
        MatrixXd InitialState = MatrixXd(phis.size(), 1);
        TransformationMatrix = MatrixXcd(phis.size(), phis.size());
        for (auto it = phis.begin(); it != phis.end(); ++it) {
            PHINode* PN = *it;
            try {
                InitialState(count, 0) = DoubleFromValue(PN->getIncomingValue(0));
            } catch (int) {
                return false;
            }
            if (!trackUpdates(count, PN)) {
                return false;
            }
            ++count;
        }

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
        
        /* XXX:
         * return true;
         */
        return false;
    }

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

    bool trackUpdates(int count, PHINode* PN) {
        return false;
    }

    int64_t IntFromValue(Value* V) {
        return cast<ConstantInt>(V)->getValue().getSExtValue();
    }

    double DoubleFromValue(Value* V) {
        if (isa<ConstantInt>(V)) {
            return double(IntFromValue(V));
        } else if (isa<ConstantFP>(V)) {
            return cast<ConstantFP>(V)->getValueAPF().convertToDouble();
        } else {
            throw 0;
        }
    }
};

char ADPass::ID = 0;

}
 
static RegisterPass<ADPass> X("auto-diagonalize",
    "Diagonalize linear dynamical systems", false, false);

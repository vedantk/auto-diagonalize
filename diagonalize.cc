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
    std::vector<BasicBlock*> blocks;

    /* Keep track of the loop iterator and state variables. */
    int iter_base;
    Value* iter_final;
    PHINode* iter_var;

    static char ID;

    ADPass() : LoopPass(ID) {}

    virtual void getAnalysisUsage(AnalysisUsage&) const {}

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

                errs() << *irhs->getParent() << "\n";

                /* loop->contains(irhs) causes bizarre segfaults */
                if (!(std::find(blocks.begin(), blocks.end(),
                    irhs->getParent()) == blocks.end()))
                {
                    return false;
                }
            }

            Value *incrLhs, *incrRhs;
            if (match(cmpLhs, m_Add(m_Value(incrLhs), m_Value(incrRhs)))) {
                if (isa<PHINode>(incrLhs)
                    && isa<ConstantInt>(incrRhs)
                    && cast<ConstantInt>(incrRhs)->isOne())
                {
                    iter_var = cast<PHINode>(incrLhs);
                    iter_final = cmpRhs;
                    errs() << "all good! in extract\n";
                    return true;
                }
            }
        }
        errs() << "oh shit\n";
        return false;
    }

    virtual bool runOnLoop(Loop* loop, LPPassManager&) {
        if (loop->getSubLoops().size()) {
            return false;
        }

        if (!loop->getUniqueExitBlock()) {
            return false;
        }

        blocks = loop->getBlocks();

        errs() << "***** Loop contents *****\n";
        for (size_t i=0; i < blocks.size(); ++i) {
            errs() << *blocks[i] << "\n";
        }
        errs() << "***** Loop contents *****\n";

        /* Filter away loops with invalid instructions. */
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
                    if (PN->getNumIncomingValues() != 2) {
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

        /* Extract the induction variable if possible. */
        if (!extractIterator(icmp)) {
            return false;
        } else {
            phis.erase(iter_var);
        }

        for (auto it = phis.begin(); it != phis.end(); ++it) {
            errs() << **it << "\n";
        }

        return false;
    }
};

char ADPass::ID = 0;

}
 
static RegisterPass<ADPass> X("auto-diagonalize",
    "Diagonalize linear dynamical systems", false, false);

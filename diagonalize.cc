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
	DominatorTree* DT; /* XXX */
	std::vector<BasicBlock*> blocks;

	/* Keep track of the loop iterator and state variables. */
	int iter_base;
	Value* iter_final;
	PHINode* iter_var;

	static char ID;

	ADPass() : LoopPass(ID) {}

	virtual void getAnalysisUsage(AnalysisUsage& AU) const {
		AU.addRequired<DominatorTree>();
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
			Value *incrLhs, *incrRhs;
			if (match(cmpLhs, m_Add(m_Value(incrLhs), m_Value(incrRhs)))) {
				if (isa<PHINode>(incrLhs)
					&& isa<ConstantInt>(incrRhs)
					&& cast<ConstantInt>(incrRhs)->isOne()
					&& cmpRhs->isUsedInBasicBlock(loop->getLoopPredecessor()))
				{
					iter_var = cast<PHINode>(incrLhs);
					iter_final = cmpRhs;
					return true;
				}
			}
		}
		return false;
	}

	virtual bool runOnLoop(Loop* loop, LPPassManager&) {
		errs() << "In runOnLoop\n";

		if (loop->getSubLoops().size()) {
			errs() << "subloops = " << loop->getSubLoops().size() << "\n";
			return false;
		}

		if (!loop->getUniqueExitBlock()) {
			errs() << "no unique exit block\n";
			return false;
		}

		blocks = loop->getBlocks();

		errs() << "***** Loop contents *****\n";
		for (auto it = blocks.begin(); it != blocks.end(); ++it) {
			errs() << *it << "\n";
		}
		errs() << "***** Loop contents *****\n";

		errs() << "Basic sanity checks done\n";

		/* Filter away loops with invalid instructions. */
		int nr_cmps = 0;
		ICmpInst* icmp = NULL;
		SmallPtrSet<PHINode*, 4> phis;
		for (auto BI = blocks.begin(); BI != blocks.end(); ++BI) {
			BasicBlock* BB = *BI;
			for (auto II = BB->begin(); II != BB->end(); ++II) {
				Instruction* instr = static_cast<Instruction*>(II);
				if (isa<ICmpInst>(instr)) {
					if (++nr_cmps > 1) {
						errs() << "Found too many cmps\n";
						return false;
					} else {
						icmp = cast<ICmpInst>(instr);
					}
				} else if (isa<PHINode>(instr)) {
					PHINode* PN = cast<PHINode>(instr);
					if (PN->getNumIncomingValues() != 2) {
						errs() << "PHI " << PN << " has too in edges\n";
						return false;
					} else {
						phis.insert(PN);
					}
				} else if (isa<BinaryOperator>(instr)) {
					BinaryOperator* binop = cast<BinaryOperator>(instr);
					if (!(binop->getOpcode() >= Instruction::Add
						&& binop->getOpcode() <= Instruction::FDiv))
					{
						errs() << "Blacklisted binop; " << binop << "\n";
						return false;
					}
				} else if (!(isa<BranchInst>(instr)
							|| isa<IndirectBrInst>(instr)))
				{
					errs() << "Bizarre unexpected instr; " << instr << "\n";
					return false;
				}
			}
		}

		/* Extract the induction variable if possible. */
		DT = &getAnalysis<DominatorTree>(); /* XXX */
		if (icmp == NULL || !extractIterator(icmp)) {
			return false;
		} else {
			phis.erase(iter_var);
		}

		for (auto it = phis.begin(); it != phis.end(); ++it) {
			errs() << *it << "\n";
		}

		return false;
	}
};

char ADPass::ID = 0;

}
 
static RegisterPass<ADPass> X("auto-diagonalize",
	"Diagonalize linear dynamical systems", false, false);

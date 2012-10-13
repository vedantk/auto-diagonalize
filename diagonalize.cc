/*
 * diagonalize.cc
 */

#include "llvm/Pass.h"
#include "llvm/Function.h"
#include "llvm/Instructions.h"
#include "llvm/ValueSymbolTable.h"
#include "llvm/Analysis/LoopInfo.h"

#include "llvm/ADT/StringMap.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/Support/raw_ostream.h"

using namespace llvm;

namespace
{

typedef SmallVector<BasicBlock*, 4> BlockVector;

struct ADPass : public FunctionPass {
	Loop* loop;
	LoopInfo* LI;
	BasicBlock *loop_exit, *loop_hdr;
	StringMap<bool> sysVars, loopSysVars;

	static char ID;

	ADPass() : FunctionPass(ID) {}

	virtual void getAnalysisUsage(AnalysisUsage& AU) const {
		AU.addRequired<LoopInfo>();
	}

	virtual bool runOnFunction(Function& F) {
		/* XXX */
		// errs() << F;
		/* XXX */

		bool changed = false;
		LI = &getAnalysis<LoopInfo>();

		ValueSymbolTable fvals = F.getValueSymbolTable();
		for (auto it = fvals.begin(); it != fvals.end(); ++it) {
			/* Initially, we optimistically assume that every
			 * variable in the function symtable is a parameter. */
			sysVars.GetOrCreateValue(it->getKey(), true);
		}

		for (auto bb = F.begin(), e = F.end(); bb != e; ++bb) {
			if (!(loop = LI->getLoopFor(bb))) {
				continue;
			}

			if (loop->getLoopDepth() != 1) {
				continue;
			}

			if (!LI->isLoopHeader(bb)) {
				continue;
			}
			loop_hdr = bb;

			/* Only handle loops with one exit point. */
			if (!(loop_exit = loop->getExitBlock())) {
				continue;
			}

			BlockVector lextent;
			lextent.push_back(loop_hdr);

			/* Our loop must be contiguous, and nest nothing. */
			++bb;
			bool defectiveLoop = false;
			while (static_cast<BasicBlock*>(bb) != loop_exit) {
				if (LI->getLoopFor(bb) != loop) {
					defectiveLoop = true;
					break;
				}
				lextent.push_back(bb);
				++bb;
			}
			if (defectiveLoop) {
				continue;
			}
			lextent.push_back(loop_exit);

			/* Copy the Function's sysVars: not every toplevel
			 * loop in the function must share the same state. */
			loopSysVars = StringMap<bool>(sysVars);
			changed |= processLoop(lextent);
		}

		return changed;
	}

	bool processLoop(BlockVector& lextent) {
		for (auto bb = lextent.begin(); bb != lextent.end(); ++bb) {
			BasicBlock* blk = *bb;

			/* XXX */
			errs() << *blk << "\n";
			/* XXX */

			for (auto it = blk->begin(); it != blk->end(); ++it) {
				if (sysInconsistent(it)) {
					return false;
				}
			}
		}

		errs() << "------------------------------------------\n";
		return false;
	}

	bool sysInconsistent(Instruction* instr) {
		if (isa<PHINode>(instr)) {

		} else if (isa<ReturnInst>(instr)) {

		} else if (isa<BinaryOperator>(instr)) {

		}
		return true;
	}
};

char ADPass::ID = 0;

}
 
static RegisterPass<ADPass> X("auto-diagonalize",
	"Diagonalize linear dynamical systems", false, false);

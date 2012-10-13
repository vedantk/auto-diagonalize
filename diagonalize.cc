/*
 * diagonalize.cc
 */

#include "llvm/Pass.h"
#include "llvm/Function.h"
#include "llvm/Constants.h"
#include "llvm/Instructions.h"
#include "llvm/ValueSymbolTable.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/ADT/ValueMap.h"
#include "llvm/Support/raw_ostream.h"

using namespace llvm;

namespace
{

struct ADPass : public FunctionPass {
	Loop* loop;
	LoopInfo* LI;
	BasicBlock* backedge; /* -> loop_hdr */
	BasicBlock* loop_hdr; /* -> loop_exit */
	BasicBlock* loop_exit; /* (-> loop_hdr) | (-> exit) */
	Constant* iter_base; /* Initial value of the iterator. */
	ValueMap<Value*, Constant*> x0; /* Initial system state. */

	static char ID;

	ADPass() : FunctionPass(ID) {}

	virtual void getAnalysisUsage(AnalysisUsage& AU) const {
		AU.addRequired<LoopInfo>();
	}

	virtual bool runOnFunction(Function& F) {
		LI = &getAnalysis<LoopInfo>();
		ValueSymbolTable& symtab = F.getValueSymbolTable();

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

			/* Only handle loops with a header and one backedge. */
			++bb;
			if (static_cast<BasicBlock*>(bb) != loop_exit) {
				continue;
			}

			/* XXX */
			errs() << *loop_hdr << *loop_exit << "\n";

			if (transformLoop(symtab)) {
				return true;
			}
		}

		return false;
	}

	bool transformLoop(ValueSymbolTable& symtab) {
		backedge = *pred_begin(loop_hdr);
		for (auto kv = symtab.begin(); kv != symtab.end(); ++kv) {
			Value* val = kv->getValue();

			/* XXX */
			errs() << kv->getKey() << "\n";
			errs() << ":: " << *val << "\n";

			if (isa<PHINode>(val)) {
				PHINode* PN = cast<PHINode>(val);

				/* We only want to deal with simple phi nodes. */
				if (PN->getNumIncomingValues() != 2) {
					continue;
				}

				checkIncomingBlock(PN, PN->getIncomingBlock(0));
				checkIncomingBlock(PN, PN->getIncomingBlock(1));
			}
		}

		if (checkCmps(loop_hdr) || checkCmps(loop_exit)) {
			return false;
		}

		return false;
	}

	void checkIncomingBlock(PHINode* PN, BasicBlock* incoming) {
		/* If the phi has an incoming value from outside of the loop,
		 * we can tentatively call it a system parameter. */
		if (incoming != loop_hdr && incoming != loop_exit) {
			Value* val = PN->getIncomingValueForBlock(incoming);
			if (isa<Constant>(val)) {
				if (!x0.count(PN)) {
					x0[PN] = cast<Constant>(val);
				} else {
					/* If there is a conflicting constant, give up. */
					x0.erase(PN);
				}
			}
		}
	}

	bool checkCmps(BasicBlock* blk) {
		/* System parameters may not be in comparisons, though the loop
		 * iterator should be. Check if the iterator is sane. */
		bool hasiter = false;
		for (auto it = blk->begin(); it != blk->end(); ++it) {
			Instruction* instr = static_cast<Instruction*>(it);
			if (isa<CmpInst>(instr)) {
				hasiter |= checkCmpOperand(instr->getOperand(0));
				hasiter |= checkCmpOperand(instr->getOperand(1));
			}
		}
		return hasiter;
	}

	bool checkCmpOperand(Value* val) {
		/* Determine whether or not this CMP operand is an iterator, and
		 * remove it from the initial system state if it's in there. */
		PHINode* CPN;
		if (isa<PHINode>(val)) {
			CPN = cast<PHINode>(val);
		} else {
			return false;
		}
		
		if (x0.count(CPN)) {
			iter_base = x0[CPN];
			x0.erase(CPN);
		}

		if (Instruction* Ins =
			dyn_cast<Instruction>(CPN->getIncomingValueForBlock(backedge)))
		{
			if (Ins->getOpcode() == Instruction::Add
				&& Ins->getOperand(0) == CPN)
			{
				if (ConstantInt* CI =
					dyn_cast<ConstantInt>(Ins->getOperand(1)))
				{
					if (CI->equalsInt(1)) {
						return true;
					}
				}
			}
		}
		return false;
	}
};

char ADPass::ID = 0;

}
 
static RegisterPass<ADPass> X("auto-diagonalize",
	"Diagonalize linear dynamical systems", false, false);

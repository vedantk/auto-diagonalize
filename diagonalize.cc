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
	DominatorTree* DT;
	BasicBlock* backedge; /* -> loop_hdr */
	BasicBlock* loop_hdr; /* -> loop_exit */
	BasicBlock* loop_exit; /* (-> loop_hdr) | (-> exit) */
	Constant* iter_base; /* Initial value of the iterator. */
	ValueMap<Value*, Constant*> x0; /* Initial system state. */

	static char ID;

	ADPass() : FunctionPass(ID) {}

	virtual void getAnalysisUsage(AnalysisUsage& AU) const {
		AU.addRequired<LoopInfo>();
		AU.addRequired<DominatorTree>();
	}

	virtual bool runOnFunction(Function& F) {
		errs() << F;

		x0.clear();
		LI = &getAnalysis<LoopInfo>();
		DT = &getAnalysis<DominatorTree>();
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

			if (isa<PHINode>(val)) {
				PHINode* PN = cast<PHINode>(val);
				if (PN->getParent() != loop_hdr) {
					continue;
				}

				/* We only want to deal with simple phi nodes. */
				if (PN->getNumIncomingValues() != 2) {
					continue;
				}

				checkIncomingBlock(PN, PN->getIncomingBlock(0));
				checkIncomingBlock(PN, PN->getIncomingBlock(1));
			}
		}

		if (!checkCmps(loop_hdr) && !checkCmps(loop_exit)) {
			return false;
		}

		for (auto kv = x0.begin(); kv != x0.end(); ++kv) {
			errs() << *kv->first << " : " << *kv->second << "\n";
		}

		return false;
	}

	void checkIncomingBlock(PHINode* PN, BasicBlock* incoming) {
		/* If the phi has an incoming value from outside of the loop,
		 * we can tentatively call it a system parameter. */
		if (DT->dominates(incoming, loop_hdr)) {
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

	PHINode* traverseToKnownPhi(Value* val) {
		/* Given an operand of an instruction, traverse up the basic block
		 * until you can resolve the value to a known PHI node. */
		if (val == NULL) {
			return NULL;
		} else if (isa<PHINode>(val)) {
			PHINode* VPN = cast<PHINode>(val);
			if (x0.count(VPN)) {
				return VPN;
			}
		} else if (isa<Instruction>(val)) {
			Instruction* ins = cast<Instruction>(val);

			/* Be careful not to traverse too far up. */
			if (DT->properlyDominates(ins->getParent(), loop_hdr)) {
				return NULL;
			}

			for (size_t i=0; i < ins->getNumOperands(); ++i) {
				PHINode* opval = traverseToKnownPhi(ins->getOperand(i));
				if (opval) {
					return opval;
				}
			}
		}
		return NULL;
	}

	bool checkCmpOperand(Value* val) {
		/* Determine whether or not this CMP operand is an iterator, and
		 * remove it from the initial system state if it's in there. */
		PHINode* CPN = traverseToKnownPhi(val);
		if (CPN == NULL) {
			return false;
		}

		iter_base = x0[CPN];
		x0.erase(CPN);

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

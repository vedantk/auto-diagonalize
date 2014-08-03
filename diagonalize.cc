/*
* Copyright (c) 2013 Vedant Kumar <vsk@berkeley.edu>
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software. THE SOFTWARE IS
* PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
* INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
* FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*/

/*
 * diagonalize.cc
 */

#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/Instructions.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/ADT/ValueMap.h"
#include "llvm/Support/PatternMatch.h"
#include "llvm/Support/InstIterator.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

using namespace llvm;
using namespace PatternMatch;
using namespace Eigen;

#define OP_IN_RANGE(_op, _start, _end) \
    (_op >= Instruction:: _start && _op <= Instruction:: _end)

namespace
{

class ADPass : public LoopPass
{
    Loop* loop;
    DominatorTree* DT;
    ICmpInst* loop_cond; 
    BasicBlock* loop_body;
    BasicBlock* exit_block;
    Value* nr_iters;
    BasicBlock* dgen; 
    ConstantInt* start_iter;
    PHINode* iter_var;
    ValueMap<PHINode*, size_t> phis;

    /* Track linear combinations of state variables. */
    typedef ValueMap<PHINode*, double> Coefficients;

public:
    static char ID;

    ADPass() : LoopPass(ID) {}

    virtual void getAnalysisUsage(AnalysisUsage& AU) const {
        AU.addRequired<DominatorTree>();
    }

    /*
     * Drive the optimization pass.
     * Return true iff the optimization was applied.
     */
    virtual bool runOnLoop(Loop* L, LPPassManager&) {
        loop = L;
        DT = &getAnalysis<DominatorTree>();
        phis.clear();

        if (loop->getBlocks().size() != 1
            || !(loop_body = loop->getBlocks().front())
            || loop->getSubLoops().size()
            || !(exit_block = loop->getUniqueExitBlock())
            || !loopFilter()
            || !extractIterator())
        {
            return false;
        }

        /* Assign a dimension to each state variable. */
        size_t phi_index = 0;
        phis.erase(iter_var);
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
            InitialState(phi_label, 0) = toDouble(getPhiConstVal(PN));

            Coefficients coeffs;
            if (!trackUpdates(getPhiFeedbackVal(PN), coeffs)) {
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
        LLVMContext& ctx = loop_body->getContext();
        Function* parentFunc = exit_block->getParent();

        Module* mod = parentFunc->getParent();
        dgen = BasicBlock::Create(ctx, "dgen", parentFunc, exit_block);
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
        Value* iexpt = BinaryOperator::Create(Instruction::Sub, nr_iters,
           start_iter, "iexpt", dgen);
        iexpt = BinaryOperator::Create(Instruction::Add, iexpt,
           ConstantInt::get(Type::getInt32Ty(ctx), 1), "iexpt_adj", dgen);
        Value* exponent = CastInst::Create(Instruction::UIToFP, iexpt, numTy,
            "fexpt", dgen);
        for (size_t j = 0; j < nr_phis; ++j) {
            /* λ[j]^n */
            Value* exptargs[] = {
                toConstantFP(ctx, std::real(D(j, j))), exponent
            };
            Value* eigvexpt = CallInst::Create(powf,
                ArrayRef<Value*>(exptargs, 2), "eigvexpt", dgen);

            /* PDn[i][j] = P[i][j] * λ[j]^n */
            for (size_t i = 0; i < nr_phis; ++i) {
                size_t index = i * nr_phis + j;
                PDn[index] = BinaryOperator::Create(Instruction::FMul,
                    toConstantFP(ctx, std::real(P(i, j))), eigvexpt, "pdn",
                                 dgen);
            }
        }

        /* xf = P(D^n) * Pinv * x0 */
        Value** soln = new Value*[nr_phis];
        Value* zero = toConstantFP(ctx, 0.0);
        for (size_t i = 0; i < nr_phis; ++i) {
            soln[i] = zero;
            for (size_t j = 0; j < nr_phis; ++j) {
                /* dotp = <a, b> : a ∈ row(P(D^n)), b ∈ col(Pinv) */
                Value* dotp = zero;
                for (size_t k = 0; k < nr_phis; ++k) {
                    Value* ik_kj = BinaryOperator::Create(Instruction::FMul,
                        PDn[i * nr_phis + k],
                        toConstantFP(ctx, std::real(Pinv(k, j))),
                        "ik_kj", dgen);
                    dotp = BinaryOperator::Create(Instruction::FAdd,
                        ik_kj, dotp, "dotp", dgen);
                }

                /* xf[i] = ∑ P(D^n)Pinv[i][j] * x0[j] */
                Value* xj_prod = BinaryOperator::Create(Instruction::FMul,
                    dotp, toConstantFP(ctx, InitialState(j)), "pdpxj",
                    dgen);
                soln[i] = BinaryOperator::Create(Instruction::FAdd, 
                    xj_prod, soln[i], "xf", dgen);
            }
        }

        delete[] PDn;

        /* Rewire edges headed in and out of the loop. */
        for (inst_iterator II = inst_begin(parentFunc),
                           E = inst_end(parentFunc); II != E; ++II)
        {
            Instruction* instr = &*II;
            if (!isa<BranchInst>(instr)) {
                continue;
            }

            BranchInst* BI = cast<BranchInst>(instr);
            for (unsigned k = 0; k < BI->getNumSuccessors(); ++k) {
                if (BI->getSuccessor(k) == loop_body) {
                    BI->setSuccessor(k, dgen);
                }
            }
        }
        BranchInst::Create(exit_block, dgen);

        /* Replace values leading into phi nodes. */
        for (auto kv = phis.begin(); kv != phis.end(); ++kv) {
            PHINode* loopPhi = kv->first;
            Value* incoming = getPhiFeedbackVal(loopPhi);
            Value* target = soln[kv->second];
            rewriteLiveValues(incoming, target);
        }
        loop_body->replaceSuccessorsPhiUsesWith(dgen);

        /* Delete the old loop. */
        loop_body->eraseFromParent();

        delete[] soln;

        return true;
    }

private:
    /*
     * Check whether the loop is linearizale.
     */
    bool loopFilter() {
        loop_cond = NULL;
        for (auto II = loop_body->begin(); II != loop_body->end(); ++II) {
            Instruction* instr = II;
            if (isa<ICmpInst>(instr)) {
                if (loop_cond) {
                    return false;
                } else {
                    loop_cond = cast<ICmpInst>(instr);
                }
            } else if (isa<PHINode>(instr)) {
                PHINode* PN = cast<PHINode>(instr);
                if (PN->getNumIncomingValues() != 2
                    || !DT->dominates(PN->getIncomingBlock(0), loop_body))
                {
                    return false;
                }
                Value* inLhs = PN->getIncomingValue(0);
                Value* inRhs = PN->getIncomingValue(1);
                if (!isa<ConstantInt>(inLhs) && !isConstant(inRhs)) {
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
        return loop_cond != NULL;
    }

    /*
     * Given a loop in canonical form, extract the loop condition and the
     * starting iteration. Return true iff basic sanity checks pass.
     */
    bool extractIterator() {
        Value* loop_var = NULL;
        ICmpInst::Predicate IPred;
        if (!match(loop_cond, m_ICmp(IPred, m_Value(loop_var),
                                            m_Value(nr_iters)))
            || IPred != CmpInst::Predicate::ICMP_EQ
            || !isa<PHINode>(loop_var))
        {
            return false;
        }

        iter_var = cast<PHINode>(loop_var);
        if (!match(getPhiConstVal(iter_var), m_ConstantInt(start_iter))) {
            return false;
        }

        bool foundIncr = false;
        for (auto II = loop_body->begin(); II != loop_body->end(); ++II) {
            Instruction* instr = II;
            if (!isa<BinaryOperator>(II)) {
                continue;
            }

            BinaryOperator* binop = cast<BinaryOperator>(instr);
            if (foundIncr &&
                (binop->getOperand(0) == iter_var
                || binop->getOperand(1) == iter_var))
            {
                return false;
            }

            if (binop->getOpcode() == Instruction::Add
                && ((binop->getOperand(0) == iter_var
                      && isa<ConstantInt>(binop->getOperand(1))
                      && toInt(binop->getOperand(1)) == 1) ||
                    (binop->getOperand(1) == iter_var
                      && isa<ConstantInt>(binop->getOperand(0))
                      && toInt(binop->getOperand(0)) == 1)))
            {
                foundIncr = true;
            }
        }
        return foundIncr;
    }

    Value* getPhiConstVal(PHINode* PN) {
      if (isConstant(PN->getIncomingValue(0))) {
        return PN->getIncomingValue(0);
      }
      return PN->getIncomingValue(1);
    }

    Value* getPhiFeedbackVal(PHINode* PN) {
      if (isConstant(PN->getIncomingValue(0))) {
        return PN->getIncomingValue(1);
      }
      return PN->getIncomingValue(0);
    }

    /*
     * Find a linear combination of state variables which generate @parent.
     * Return true iff a valid list of coefficients is found.
     */
    bool trackUpdates(Value* parent, Coefficients& coeffs, bool root = true) {
        if (isa<Constant>(parent)) {
            /*
             * If a state variable is set to a constant during each iteration,
             * the compiler should lift it out of the loop before we get here.
             */
            return !root;
        } else if (isa<PHINode>(parent)) {
            PHINode* PN = cast<PHINode>(parent);
            coeffs[PN] = 1;
            return phis.count(PN) == 1;
        } else if (isa<BinaryOperator>(parent)) {
            BinaryOperator* binop = cast<BinaryOperator>(parent);
            int opcode = binop->getOpcode();
            Value *LHS = binop->getOperand(0),
                  *RHS = binop->getOperand(1);

            Coefficients lhsCoeffs, rhsCoeffs;
            if (!trackUpdates(LHS, lhsCoeffs, false)
                || !trackUpdates(RHS, rhsCoeffs, false))
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
                /* Mul instructions should only have one scalar operand. */
                if (!(isScalar(lhsCoeffs) ^ isScalar(rhsCoeffs))) {
                    return false;
                }

                /* Div instructions cannot have scalar numerators. */
                if (OP_IN_RANGE(opcode, UDiv, FDiv) && isScalar(lhsCoeffs)) {
                    return false;
                }

                scalar = toDouble(isScalar(lhsCoeffs) ? LHS : RHS);
            }

            /* Merge the two sets of coefficients. */
            for (auto kv = phis.begin(); kv != phis.end(); ++kv) {
                PHINode* PN = kv->first;
                double lcoeff = lhsCoeffs.lookup(PN);
                double rcoeff = rhsCoeffs.lookup(PN);

                /* Adding nil entries to 'coeffs' breaks isScalar(). */
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
            return true;
        }
        return false;
    }

    bool isConstant(Value* V) {
      return isa<ConstantInt>(V) || isa<ConstantFP>(V);
    }

    int64_t toInt(Value* V) {
        return cast<ConstantInt>(V)->getValue().getSExtValue();
    }

    double toDouble(Value* V) {
        if (isa<ConstantInt>(V)) {
            return double(toInt(V));
        } else if (isa<ConstantFP>(V)) {
            return cast<ConstantFP>(V)->getValueAPF().convertToDouble();
        }
        return 0.0;
    }

    Value* toConstantFP(LLVMContext& ctx, double n) {
        return ConstantFP::get(ctx, APFloat(n));
    }

    bool isScalar(Coefficients& coeffs) {
        return coeffs.size() == 0;
    }

    /*
     * Verify that A == PDP^-1 within an acceptable margin of error.
     */
    bool checkSystem(MatrixXd& A, MatrixXcd& P,
                     MatrixXcd& D, MatrixXcd& Pinv)
    {
        MatrixXcd PDPi = P*D*Pinv;
        const double epsilon = 25 * std::numeric_limits<double>::epsilon();
        for (int i = 0; i < P.rows(); ++i) {
            for (int j = 0; j < P.cols(); ++j) {
                if (std::imag(P(i, j)) != 0.0) {
                    return false;
                }
                if (std::abs(A(i, j) - PDPi(i, j)) > epsilon) {
                    return false;
                }
            }
        }
        return true;
    }

    /* 
     * Some instructions may use the results of the loop. Rewire these
     * instructions so that they use the values from the @dgen basic block.
     */
    void rewriteLiveValues(Value* oldval, Value* target) {
        for (auto UI = oldval->use_begin(); UI != oldval->use_end(); ++UI) {
            User* user = *UI;

            /* Don't bother updating values in @loop_body. It dies soon. */
            if (isPartialComputation(user)) {
                continue;
            }

            if (isa<PHINode>(user)) {
                PHINode* userPhi = cast<PHINode>(user);
                for (unsigned i=0; i < userPhi->getNumIncomingValues(); ++i) {
                    if (userPhi->getIncomingValue(i) == oldval) {
                        userPhi->setIncomingValue(i, target);
                        userPhi->setIncomingBlock(i, dgen);
                    }
                }
                continue;
            }

            for (auto op = user->op_begin(); op != user->op_end(); ++op) {
                Use& phiUse = *op;
                Value* operand = phiUse.get();
                if (operand == oldval) {
                    phiUse.set(target);
                }
            }
        }
    }

    /* 
     * Check if the target value is defined in the loop body.
     */
    bool isPartialComputation(Value* target) {
        for (auto II = loop_body->begin(); II != loop_body->end(); ++II) {
            Instruction* instr = II;
            if (instr == target) {
                return true;
            }
        }
        return false;
    }
};

char ADPass::ID = 0;

}

static RegisterPass<ADPass> X("auto-diagonalize",
    "Diagonalize linear dynamical systems", false, false);

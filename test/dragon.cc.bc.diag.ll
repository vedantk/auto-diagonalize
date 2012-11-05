; ModuleID = 'test/dragon.cc.bc.diag'
target datalayout = "e-p:64:64:64-i1:8:8-i8:8:8-i16:16:16-i32:32:32-i64:64:64-f32:32:32-f64:64:64-v64:64:64-v128:128:128-a0:0:64-s0:64:64-f80:128:128-n8:16:32:64-S128"
target triple = "x86_64-unknown-linux-gnu"

@.str = private unnamed_addr constant [25 x i8] c"Usage: ./dragon <number>\00", align 1
@.str1 = private unnamed_addr constant [17 x i8] c"dragon(%d) = %f\0A\00", align 1

define double @_Z6dragoni(i32 %n) nounwind uwtable readnone {
  %1 = icmp slt i32 %n, 3
  br i1 %1, label %._crit_edge, label %.lr.ph

.lr.ph:                                           ; preds = %0
  %2 = add i32 %n, 1
  br label %dgen

dgen:                                             ; preds = %.lr.ph
  %iexpt = sub i32 %2, 3
  %fexpt = uitofp i32 %iexpt to double
  %eigvexpt = call double @llvm.pow.f64(double 0xBFBB0404CAEA1034, double %fexpt)
  %pdn = fmul double 0xBFEC1D6FA818187F, %eigvexpt
  %pdn1 = fmul double 0x3FD96E64D270F57C, %eigvexpt
  %pdn2 = fmul double 0x3FD0F44336F5F8FF, %eigvexpt
  %eigvexpt3 = call double @llvm.pow.f64(double 0x3FF1B0404CAEA101, double %fexpt)
  %pdn4 = fmul double 0x3FC623BE35035B23, %eigvexpt3
  %pdn5 = fmul double 0x3FEA396BF8F24F62, %eigvexpt3
  %pdn6 = fmul double 0x3FE17B9D50A18A42, %eigvexpt3
  %eigvexpt7 = call double @llvm.pow.f64(double 1.000000e+00, double %fexpt)
  %pdn8 = fmul double 0xBCC23BACD5746AD9, %eigvexpt7
  %pdn9 = fmul double 0xBFEC9F25C5BFEDE1, %eigvexpt7
  %pdn10 = fmul double 0x3FDC9F25C5BFEDB4, %eigvexpt7
  %ik_kj = fmul double %pdn, 0xBFF09FC209315135
  %dotp = fadd double %ik_kj, 0.000000e+00
  %ik_kj11 = fmul double %pdn4, 0x3FE01F0E0F51E380
  %dotp12 = fadd double %ik_kj11, %dotp
  %ik_kj13 = fmul double %pdn8, 0x3C8EAA83E629112D
  %dotp14 = fadd double %ik_kj13, %dotp12
  %pdpxj = fmul double %dotp14, 1.000000e+00
  %xf = fadd double %pdpxj, 0.000000e+00
  %ik_kj15 = fmul double %pdn, 0x3FB80F3A96F66DB6
  %dotp16 = fadd double %ik_kj15, 0.000000e+00
  %ik_kj17 = fmul double %pdn4, 0x3FDE8D9247933F9A
  %dotp18 = fadd double %ik_kj17, %dotp16
  %ik_kj19 = fmul double %pdn8, 0xBFE471AD441B60C9
  %dotp20 = fadd double %ik_kj19, %dotp18
  %pdpxj21 = fmul double %dotp20, 2.000000e+00
  %xf22 = fadd double %pdpxj21, %xf
  %ik_kj23 = fmul double %pdn, 0x3FC80F3A96F66DB5
  %dotp24 = fadd double %ik_kj23, 0.000000e+00
  %ik_kj25 = fmul double %pdn4, 0x3FEE8D9247933FD2
  %dotp26 = fadd double %ik_kj25, %dotp24
  %ik_kj27 = fmul double %pdn8, 0x3FEEAA83E629112D
  %dotp28 = fadd double %ik_kj27, %dotp26
  %pdpxj29 = fmul double %dotp28, 4.000000e+00
  %xf30 = fadd double %pdpxj29, %xf22
  %ik_kj31 = fmul double %pdn1, 0xBFF09FC209315135
  %dotp32 = fadd double %ik_kj31, 0.000000e+00
  %ik_kj33 = fmul double %pdn5, 0x3FE01F0E0F51E380
  %dotp34 = fadd double %ik_kj33, %dotp32
  %ik_kj35 = fmul double %pdn9, 0x3C8EAA83E629112D
  %dotp36 = fadd double %ik_kj35, %dotp34
  %pdpxj37 = fmul double %dotp36, 1.000000e+00
  %xf38 = fadd double %pdpxj37, 0.000000e+00
  %ik_kj39 = fmul double %pdn1, 0x3FB80F3A96F66DB6
  %dotp40 = fadd double %ik_kj39, 0.000000e+00
  %ik_kj41 = fmul double %pdn5, 0x3FDE8D9247933F9A
  %dotp42 = fadd double %ik_kj41, %dotp40
  %ik_kj43 = fmul double %pdn9, 0xBFE471AD441B60C9
  %dotp44 = fadd double %ik_kj43, %dotp42
  %pdpxj45 = fmul double %dotp44, 2.000000e+00
  %xf46 = fadd double %pdpxj45, %xf38
  %ik_kj47 = fmul double %pdn1, 0x3FC80F3A96F66DB5
  %dotp48 = fadd double %ik_kj47, 0.000000e+00
  %ik_kj49 = fmul double %pdn5, 0x3FEE8D9247933FD2
  %dotp50 = fadd double %ik_kj49, %dotp48
  %ik_kj51 = fmul double %pdn9, 0x3FEEAA83E629112D
  %dotp52 = fadd double %ik_kj51, %dotp50
  %pdpxj53 = fmul double %dotp52, 4.000000e+00
  %xf54 = fadd double %pdpxj53, %xf46
  %ik_kj55 = fmul double %pdn2, 0xBFF09FC209315135
  %dotp56 = fadd double %ik_kj55, 0.000000e+00
  %ik_kj57 = fmul double %pdn6, 0x3FE01F0E0F51E380
  %dotp58 = fadd double %ik_kj57, %dotp56
  %ik_kj59 = fmul double %pdn10, 0x3C8EAA83E629112D
  %dotp60 = fadd double %ik_kj59, %dotp58
  %pdpxj61 = fmul double %dotp60, 1.000000e+00
  %xf62 = fadd double %pdpxj61, 0.000000e+00
  %ik_kj63 = fmul double %pdn2, 0x3FB80F3A96F66DB6
  %dotp64 = fadd double %ik_kj63, 0.000000e+00
  %ik_kj65 = fmul double %pdn6, 0x3FDE8D9247933F9A
  %dotp66 = fadd double %ik_kj65, %dotp64
  %ik_kj67 = fmul double %pdn10, 0xBFE471AD441B60C9
  %dotp68 = fadd double %ik_kj67, %dotp66
  %pdpxj69 = fmul double %dotp68, 2.000000e+00
  %xf70 = fadd double %pdpxj69, %xf62
  %ik_kj71 = fmul double %pdn2, 0x3FC80F3A96F66DB5
  %dotp72 = fadd double %ik_kj71, 0.000000e+00
  %ik_kj73 = fmul double %pdn6, 0x3FEE8D9247933FD2
  %dotp74 = fadd double %ik_kj73, %dotp72
  %ik_kj75 = fmul double %pdn10, 0x3FEEAA83E629112D
  %dotp76 = fadd double %ik_kj75, %dotp74
  %pdpxj77 = fmul double %dotp76, 4.000000e+00
  %xf78 = fadd double %pdpxj77, %xf70
  br label %._crit_edge

._crit_edge:                                      ; preds = %dgen, %0
  %c.0.lcssa = phi double [ 4.000000e+00, %0 ], [ %xf78, %dgen ]
  ret double %c.0.lcssa
}

define i32 @main(i32 %argc, i8** nocapture %argv) nounwind uwtable {
  %1 = icmp eq i32 %argc, 2
  br i1 %1, label %4, label %2

; <label>:2                                       ; preds = %0
  %3 = tail call i32 @puts(i8* getelementptr inbounds ([25 x i8]* @.str, i64 0, i64 0))
  br label %11

; <label>:4                                       ; preds = %0
  %5 = getelementptr inbounds i8** %argv, i64 1
  %6 = load i8** %5, align 8, !tbaa !0
  %7 = tail call i32 @atoi(i8* %6) nounwind readonly
  %8 = icmp slt i32 %7, 3
  br i1 %8, label %_Z6dragoni.exit, label %.lr.ph.i

.lr.ph.i:                                         ; preds = %4
  %9 = add i32 %7, 1
  br label %dgen

dgen:                                             ; preds = %.lr.ph.i
  %iexpt = sub i32 %9, 3
  %fexpt = uitofp i32 %iexpt to double
  %eigvexpt = call double @llvm.pow.f64(double 0xBFBB0404CAEA1034, double %fexpt)
  %pdn = fmul double 0xBFEC1D6FA818187F, %eigvexpt
  %pdn1 = fmul double 0x3FD96E64D270F57C, %eigvexpt
  %pdn2 = fmul double 0x3FD0F44336F5F8FF, %eigvexpt
  %eigvexpt3 = call double @llvm.pow.f64(double 0x3FF1B0404CAEA101, double %fexpt)
  %pdn4 = fmul double 0x3FC623BE35035B23, %eigvexpt3
  %pdn5 = fmul double 0x3FEA396BF8F24F62, %eigvexpt3
  %pdn6 = fmul double 0x3FE17B9D50A18A42, %eigvexpt3
  %eigvexpt7 = call double @llvm.pow.f64(double 1.000000e+00, double %fexpt)
  %pdn8 = fmul double 0xBCC23BACD5746AD9, %eigvexpt7
  %pdn9 = fmul double 0xBFEC9F25C5BFEDE1, %eigvexpt7
  %pdn10 = fmul double 0x3FDC9F25C5BFEDB4, %eigvexpt7
  %ik_kj = fmul double %pdn, 0xBFF09FC209315135
  %dotp = fadd double %ik_kj, 0.000000e+00
  %ik_kj11 = fmul double %pdn4, 0x3FE01F0E0F51E380
  %dotp12 = fadd double %ik_kj11, %dotp
  %ik_kj13 = fmul double %pdn8, 0x3C8EAA83E629112D
  %dotp14 = fadd double %ik_kj13, %dotp12
  %pdpxj = fmul double %dotp14, 1.000000e+00
  %xf = fadd double %pdpxj, 0.000000e+00
  %ik_kj15 = fmul double %pdn, 0x3FB80F3A96F66DB6
  %dotp16 = fadd double %ik_kj15, 0.000000e+00
  %ik_kj17 = fmul double %pdn4, 0x3FDE8D9247933F9A
  %dotp18 = fadd double %ik_kj17, %dotp16
  %ik_kj19 = fmul double %pdn8, 0xBFE471AD441B60C9
  %dotp20 = fadd double %ik_kj19, %dotp18
  %pdpxj21 = fmul double %dotp20, 2.000000e+00
  %xf22 = fadd double %pdpxj21, %xf
  %ik_kj23 = fmul double %pdn, 0x3FC80F3A96F66DB5
  %dotp24 = fadd double %ik_kj23, 0.000000e+00
  %ik_kj25 = fmul double %pdn4, 0x3FEE8D9247933FD2
  %dotp26 = fadd double %ik_kj25, %dotp24
  %ik_kj27 = fmul double %pdn8, 0x3FEEAA83E629112D
  %dotp28 = fadd double %ik_kj27, %dotp26
  %pdpxj29 = fmul double %dotp28, 4.000000e+00
  %xf30 = fadd double %pdpxj29, %xf22
  %ik_kj31 = fmul double %pdn1, 0xBFF09FC209315135
  %dotp32 = fadd double %ik_kj31, 0.000000e+00
  %ik_kj33 = fmul double %pdn5, 0x3FE01F0E0F51E380
  %dotp34 = fadd double %ik_kj33, %dotp32
  %ik_kj35 = fmul double %pdn9, 0x3C8EAA83E629112D
  %dotp36 = fadd double %ik_kj35, %dotp34
  %pdpxj37 = fmul double %dotp36, 1.000000e+00
  %xf38 = fadd double %pdpxj37, 0.000000e+00
  %ik_kj39 = fmul double %pdn1, 0x3FB80F3A96F66DB6
  %dotp40 = fadd double %ik_kj39, 0.000000e+00
  %ik_kj41 = fmul double %pdn5, 0x3FDE8D9247933F9A
  %dotp42 = fadd double %ik_kj41, %dotp40
  %ik_kj43 = fmul double %pdn9, 0xBFE471AD441B60C9
  %dotp44 = fadd double %ik_kj43, %dotp42
  %pdpxj45 = fmul double %dotp44, 2.000000e+00
  %xf46 = fadd double %pdpxj45, %xf38
  %ik_kj47 = fmul double %pdn1, 0x3FC80F3A96F66DB5
  %dotp48 = fadd double %ik_kj47, 0.000000e+00
  %ik_kj49 = fmul double %pdn5, 0x3FEE8D9247933FD2
  %dotp50 = fadd double %ik_kj49, %dotp48
  %ik_kj51 = fmul double %pdn9, 0x3FEEAA83E629112D
  %dotp52 = fadd double %ik_kj51, %dotp50
  %pdpxj53 = fmul double %dotp52, 4.000000e+00
  %xf54 = fadd double %pdpxj53, %xf46
  %ik_kj55 = fmul double %pdn2, 0xBFF09FC209315135
  %dotp56 = fadd double %ik_kj55, 0.000000e+00
  %ik_kj57 = fmul double %pdn6, 0x3FE01F0E0F51E380
  %dotp58 = fadd double %ik_kj57, %dotp56
  %ik_kj59 = fmul double %pdn10, 0x3C8EAA83E629112D
  %dotp60 = fadd double %ik_kj59, %dotp58
  %pdpxj61 = fmul double %dotp60, 1.000000e+00
  %xf62 = fadd double %pdpxj61, 0.000000e+00
  %ik_kj63 = fmul double %pdn2, 0x3FB80F3A96F66DB6
  %dotp64 = fadd double %ik_kj63, 0.000000e+00
  %ik_kj65 = fmul double %pdn6, 0x3FDE8D9247933F9A
  %dotp66 = fadd double %ik_kj65, %dotp64
  %ik_kj67 = fmul double %pdn10, 0xBFE471AD441B60C9
  %dotp68 = fadd double %ik_kj67, %dotp66
  %pdpxj69 = fmul double %dotp68, 2.000000e+00
  %xf70 = fadd double %pdpxj69, %xf62
  %ik_kj71 = fmul double %pdn2, 0x3FC80F3A96F66DB5
  %dotp72 = fadd double %ik_kj71, 0.000000e+00
  %ik_kj73 = fmul double %pdn6, 0x3FEE8D9247933FD2
  %dotp74 = fadd double %ik_kj73, %dotp72
  %ik_kj75 = fmul double %pdn10, 0x3FEEAA83E629112D
  %dotp76 = fadd double %ik_kj75, %dotp74
  %pdpxj77 = fmul double %dotp76, 4.000000e+00
  %xf78 = fadd double %pdpxj77, %xf70
  br label %_Z6dragoni.exit

_Z6dragoni.exit:                                  ; preds = %dgen, %4
  %c.0.lcssa.i = phi double [ 4.000000e+00, %4 ], [ %xf78, %dgen ]
  %10 = tail call i32 (i8*, ...)* @printf(i8* getelementptr inbounds ([17 x i8]* @.str1, i64 0, i64 0), i32 %7, double %c.0.lcssa.i)
  br label %11

; <label>:11                                      ; preds = %_Z6dragoni.exit, %2
  %.0 = phi i32 [ 1, %2 ], [ 0, %_Z6dragoni.exit ]
  ret i32 %.0
}

declare i32 @puts(i8* nocapture) nounwind

declare i32 @atoi(i8* nocapture) nounwind readonly

declare i32 @printf(i8* nocapture, ...) nounwind

declare double @llvm.pow.f64(double, double) nounwind readonly

!0 = metadata !{metadata !"any pointer", metadata !1}
!1 = metadata !{metadata !"omnipotent char", metadata !2}
!2 = metadata !{metadata !"Simple C/C++ TBAA"}

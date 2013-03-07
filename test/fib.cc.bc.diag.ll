; ModuleID = 'test/fib.cc.bc.diag'
target datalayout = "e-p:64:64:64-i1:8:8-i8:8:8-i16:16:16-i32:32:32-i64:64:64-f32:32:32-f64:64:64-v64:64:64-v128:128:128-a0:0:64-s0:64:64-f80:128:128-n8:16:32:64-S128"
target triple = "x86_64-unknown-linux-gnu"

@.str = private unnamed_addr constant [22 x i8] c"Usage: ./fib <number>\00", align 1
@.str1 = private unnamed_addr constant [15 x i8] c"fib(%d) = %lf\0A\00", align 1

define double @_Z3fibi(i32 %n) nounwind uwtable readnone {
  %1 = icmp slt i32 %n, 2
  br i1 %1, label %._crit_edge, label %.lr.ph

.lr.ph:                                           ; preds = %0
  %2 = add i32 %n, 1
  br label %dgen

dgen:                                             ; preds = %.lr.ph
  %iexpt = sub i32 %2, 2
  %fexpt = uitofp i32 %iexpt to double
  %eigvexpt = call double @llvm.pow.f64(double 0xBFE3C6EF372FE94E, double %fexpt)
  %pdn = fmul double 0xBFEB38880B4603E6, %eigvexpt
  %pdn1 = fmul double 0x3FE0D2CA0DA1530D, %eigvexpt
  %eigvexpt2 = call double @llvm.pow.f64(double 0x3FF9E3779B97F4A6, double %fexpt)
  %pdn3 = fmul double 0xBFE0D2CA0DA1530E, %eigvexpt2
  %pdn4 = fmul double 0xBFEB38880B4603E6, %eigvexpt2
  %ik_kj = fmul double %pdn, 0xBFEB38880B4603E3
  %dotp = fadd double %ik_kj, 0.000000e+00
  %ik_kj5 = fmul double %pdn3, 0xBFE0D2CA0DA1530C
  %dotp6 = fadd double %ik_kj5, %dotp
  %pdpxj = fmul double %dotp6, 1.000000e+00
  %xf = fadd double %pdpxj, 0.000000e+00
  %ik_kj7 = fmul double %pdn, 0x3FE0D2CA0DA1530D
  %dotp8 = fadd double %ik_kj7, 0.000000e+00
  %ik_kj9 = fmul double %pdn3, 0xBFEB38880B4603E5
  %dotp10 = fadd double %ik_kj9, %dotp8
  %pdpxj11 = fmul double %dotp10, 1.000000e+00
  %xf12 = fadd double %pdpxj11, %xf
  %ik_kj13 = fmul double %pdn1, 0xBFEB38880B4603E3
  %dotp14 = fadd double %ik_kj13, 0.000000e+00
  %ik_kj15 = fmul double %pdn4, 0xBFE0D2CA0DA1530C
  %dotp16 = fadd double %ik_kj15, %dotp14
  %pdpxj17 = fmul double %dotp16, 1.000000e+00
  %xf18 = fadd double %pdpxj17, 0.000000e+00
  %ik_kj19 = fmul double %pdn1, 0x3FE0D2CA0DA1530D
  %dotp20 = fadd double %ik_kj19, 0.000000e+00
  %ik_kj21 = fmul double %pdn4, 0xBFEB38880B4603E5
  %dotp22 = fadd double %ik_kj21, %dotp20
  %pdpxj23 = fmul double %dotp22, 1.000000e+00
  %xf24 = fadd double %pdpxj23, %xf18
  br label %._crit_edge

._crit_edge:                                      ; preds = %dgen, %0
  %b.0.lcssa = phi double [ 1.000000e+00, %0 ], [ %xf24, %dgen ]
  ret double %b.0.lcssa
}

define i32 @main(i32 %argc, i8** nocapture %argv) nounwind uwtable {
  %1 = icmp eq i32 %argc, 2
  br i1 %1, label %4, label %2

; <label>:2                                       ; preds = %0
  %3 = tail call i32 @puts(i8* getelementptr inbounds ([22 x i8]* @.str, i64 0, i64 0))
  br label %.loopexit

; <label>:4                                       ; preds = %0
  %5 = getelementptr inbounds i8** %argv, i64 1
  %6 = load i8** %5, align 8, !tbaa !0
  %7 = tail call i64 @strtol(i8* nocapture %6, i8** null, i32 10) nounwind
  %8 = trunc i64 %7 to i32
  %9 = icmp slt i32 %8, 2
  %10 = add i32 %8, 1
  br i1 %9, label %_Z3fibi.exit.us, label %.lr.ph.i

_Z3fibi.exit.us:                                  ; preds = %_Z3fibi.exit.us, %4
  %i.05.us = phi i32 [ %12, %_Z3fibi.exit.us ], [ 0, %4 ]
  %11 = tail call i32 (i8*, ...)* @printf(i8* getelementptr inbounds ([15 x i8]* @.str1, i64 0, i64 0), i32 %8, double 1.000000e+00)
  %12 = add nsw i32 %i.05.us, 1
  %exitcond6 = icmp eq i32 %12, 10001
  br i1 %exitcond6, label %.loopexit, label %_Z3fibi.exit.us

.lr.ph.i:                                         ; preds = %_Z3fibi.exit, %4
  %i.05 = phi i32 [ %14, %_Z3fibi.exit ], [ 0, %4 ]
  br label %dgen

dgen:                                             ; preds = %.lr.ph.i
  %iexpt = sub i32 %10, 2
  %fexpt = uitofp i32 %iexpt to double
  %eigvexpt = call double @llvm.pow.f64(double 0x3FF9E3779B97F4A6, double %fexpt)
  %pdn = fmul double 0x3FEB38880B4603E6, %eigvexpt
  %pdn1 = fmul double 0x3FE0D2CA0DA1530D, %eigvexpt
  %eigvexpt2 = call double @llvm.pow.f64(double 0xBFE3C6EF372FE94E, double %fexpt)
  %pdn3 = fmul double 0xBFE0D2CA0DA1530D, %eigvexpt2
  %pdn4 = fmul double 0x3FEB38880B4603E6, %eigvexpt2
  %ik_kj = fmul double %pdn, 0x3FEB38880B4603E4
  %dotp = fadd double %ik_kj, 0.000000e+00
  %ik_kj5 = fmul double %pdn3, 0xBFE0D2CA0DA1530C
  %dotp6 = fadd double %ik_kj5, %dotp
  %pdpxj = fmul double %dotp6, 1.000000e+00
  %xf = fadd double %pdpxj, 0.000000e+00
  %ik_kj7 = fmul double %pdn, 0x3FE0D2CA0DA1530D
  %dotp8 = fadd double %ik_kj7, 0.000000e+00
  %ik_kj9 = fmul double %pdn3, 0x3FEB38880B4603E5
  %dotp10 = fadd double %ik_kj9, %dotp8
  %pdpxj11 = fmul double %dotp10, 1.000000e+00
  %xf12 = fadd double %pdpxj11, %xf
  %ik_kj13 = fmul double %pdn1, 0x3FEB38880B4603E4
  %dotp14 = fadd double %ik_kj13, 0.000000e+00
  %ik_kj15 = fmul double %pdn4, 0xBFE0D2CA0DA1530C
  %dotp16 = fadd double %ik_kj15, %dotp14
  %pdpxj17 = fmul double %dotp16, 1.000000e+00
  %xf18 = fadd double %pdpxj17, 0.000000e+00
  %ik_kj19 = fmul double %pdn1, 0x3FE0D2CA0DA1530D
  %dotp20 = fadd double %ik_kj19, 0.000000e+00
  %ik_kj21 = fmul double %pdn4, 0x3FEB38880B4603E5
  %dotp22 = fadd double %ik_kj21, %dotp20
  %pdpxj23 = fmul double %dotp22, 1.000000e+00
  %xf24 = fadd double %pdpxj23, %xf18
  br label %_Z3fibi.exit

_Z3fibi.exit:                                     ; preds = %dgen
  %13 = tail call i32 (i8*, ...)* @printf(i8* getelementptr inbounds ([15 x i8]* @.str1, i64 0, i64 0), i32 %8, double %xf12)
  %14 = add nsw i32 %i.05, 1
  %exitcond = icmp eq i32 %14, 10001
  br i1 %exitcond, label %.loopexit, label %.lr.ph.i

.loopexit:                                        ; preds = %_Z3fibi.exit, %_Z3fibi.exit.us, %2
  %.0 = phi i32 [ 1, %2 ], [ 0, %_Z3fibi.exit ], [ 0, %_Z3fibi.exit.us ]
  ret i32 %.0
}

declare i32 @puts(i8* nocapture) nounwind

declare i32 @printf(i8* nocapture, ...) nounwind

declare i64 @strtol(i8*, i8** nocapture, i32) nounwind

declare double @llvm.pow.f64(double, double) nounwind readonly

!0 = metadata !{metadata !"any pointer", metadata !1}
!1 = metadata !{metadata !"omnipotent char", metadata !2}
!2 = metadata !{metadata !"Simple C/C++ TBAA"}

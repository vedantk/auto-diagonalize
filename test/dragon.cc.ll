; ModuleID = 'dragon.cc.bc'
target datalayout = "e-p:64:64:64-i1:8:8-i8:8:8-i16:16:16-i32:32:32-i64:64:64-f32:32:32-f64:64:64-v64:64:64-v128:128:128-a0:0:64-s0:64:64-f80:128:128-n8:16:32:64-S128"
target triple = "x86_64-unknown-linux-gnu"

@.str = private unnamed_addr constant [25 x i8] c"Usage: ./dragon <number>\00", align 1
@.str1 = private unnamed_addr constant [17 x i8] c"dragon(%d) = %f\0A\00", align 1

define double @_Z6dragoni(i32 %n) nounwind uwtable readnone {
  %1 = icmp slt i32 %n, 3
  br i1 %1, label %._crit_edge, label %.lr.ph

.lr.ph:                                           ; preds = %0
  %2 = add i32 %n, 1
  br label %3

; <label>:3                                       ; preds = %3, %.lr.ph
  %i.04 = phi i32 [ 3, %.lr.ph ], [ %11, %3 ]
  %a.03 = phi double [ 1.000000e+00, %.lr.ph ], [ %6, %3 ]
  %b.02 = phi double [ 2.000000e+00, %.lr.ph ], [ %8, %3 ]
  %c.01 = phi double [ 4.000000e+00, %.lr.ph ], [ %10, %3 ]
  %4 = fmul double %b.02, 1.000000e-01
  %5 = fmul double %c.01, 2.000000e-01
  %6 = fadd double %4, %5
  %7 = fmul double %a.03, 5.000000e-01
  %8 = fadd double %7, %b.02
  %9 = fdiv double %a.03, 3.000000e+00
  %10 = fadd double %9, %c.01
  %11 = add nsw i32 %i.04, 1
  %exitcond = icmp eq i32 %11, %2
  br i1 %exitcond, label %._crit_edge, label %3

._crit_edge:                                      ; preds = %3, %0
  %c.0.lcssa = phi double [ 4.000000e+00, %0 ], [ %10, %3 ]
  ret double %c.0.lcssa
}

define i32 @main(i32 %argc, i8** nocapture %argv) nounwind uwtable {
  %1 = icmp eq i32 %argc, 2
  br i1 %1, label %4, label %2

; <label>:2                                       ; preds = %0
  %3 = tail call i32 @puts(i8* getelementptr inbounds ([25 x i8]* @.str, i64 0, i64 0))
  br label %20

; <label>:4                                       ; preds = %0
  %5 = getelementptr inbounds i8** %argv, i64 1
  %6 = load i8** %5, align 8, !tbaa !0
  %7 = tail call i32 @atoi(i8* %6) nounwind readonly
  %8 = icmp slt i32 %7, 3
  br i1 %8, label %_Z6dragoni.exit, label %.lr.ph.i

.lr.ph.i:                                         ; preds = %4
  %9 = add i32 %7, 1
  br label %10

; <label>:10                                      ; preds = %10, %.lr.ph.i
  %i.04.i = phi i32 [ 3, %.lr.ph.i ], [ %18, %10 ]
  %a.03.i = phi double [ 1.000000e+00, %.lr.ph.i ], [ %13, %10 ]
  %b.02.i = phi double [ 2.000000e+00, %.lr.ph.i ], [ %15, %10 ]
  %c.01.i = phi double [ 4.000000e+00, %.lr.ph.i ], [ %17, %10 ]
  %11 = fmul double %b.02.i, 1.000000e-01
  %12 = fmul double %c.01.i, 2.000000e-01
  %13 = fadd double %11, %12
  %14 = fmul double %a.03.i, 5.000000e-01
  %15 = fadd double %14, %b.02.i
  %16 = fdiv double %a.03.i, 3.000000e+00
  %17 = fadd double %16, %c.01.i
  %18 = add nsw i32 %i.04.i, 1
  %exitcond.i = icmp eq i32 %18, %9
  br i1 %exitcond.i, label %_Z6dragoni.exit, label %10

_Z6dragoni.exit:                                  ; preds = %10, %4
  %c.0.lcssa.i = phi double [ 4.000000e+00, %4 ], [ %17, %10 ]
  %19 = tail call i32 (i8*, ...)* @printf(i8* getelementptr inbounds ([17 x i8]* @.str1, i64 0, i64 0), i32 %7, double %c.0.lcssa.i)
  br label %20

; <label>:20                                      ; preds = %_Z6dragoni.exit, %2
  %.0 = phi i32 [ 1, %2 ], [ 0, %_Z6dragoni.exit ]
  ret i32 %.0
}

declare i32 @puts(i8* nocapture) nounwind

declare i32 @atoi(i8* nocapture) nounwind readonly

declare i32 @printf(i8* nocapture, ...) nounwind

!0 = metadata !{metadata !"any pointer", metadata !1}
!1 = metadata !{metadata !"omnipotent char", metadata !2}
!2 = metadata !{metadata !"Simple C/C++ TBAA"}

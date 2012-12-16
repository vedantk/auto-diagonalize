; ModuleID = 'fib.cc.bc'
target datalayout = "e-p:64:64:64-i1:8:8-i8:8:8-i16:16:16-i32:32:32-i64:64:64-f32:32:32-f64:64:64-v64:64:64-v128:128:128-a0:0:64-s0:64:64-f80:128:128-n8:16:32:64-S128"
target triple = "x86_64-unknown-linux-gnu"

@.str = private unnamed_addr constant [22 x i8] c"Usage: ./fib <number>\00", align 1
@.str1 = private unnamed_addr constant [15 x i8] c"fib(%d) = %lf\0A\00", align 1

define double @_Z3fibi(i32 %n) nounwind uwtable readnone {
  %1 = icmp slt i32 %n, 2
  br i1 %1, label %._crit_edge, label %.lr.ph

.lr.ph:                                           ; preds = %0
  %2 = add i32 %n, 1
  br label %3

; <label>:3                                       ; preds = %3, %.lr.ph
  %i.03 = phi i32 [ 2, %.lr.ph ], [ %5, %3 ]
  %a.02 = phi double [ 1.000000e+00, %.lr.ph ], [ %b.01, %3 ]
  %b.01 = phi double [ 1.000000e+00, %.lr.ph ], [ %4, %3 ]
  %4 = fadd double %a.02, %b.01
  %5 = add nsw i32 %i.03, 1
  %exitcond = icmp eq i32 %5, %2
  br i1 %exitcond, label %._crit_edge, label %3

._crit_edge:                                      ; preds = %3, %0
  %b.0.lcssa = phi double [ 1.000000e+00, %0 ], [ %4, %3 ]
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
  %7 = tail call i32 @atoi(i8* %6) nounwind readonly
  %8 = icmp slt i32 %7, 2
  %9 = add i32 %7, 1
  br i1 %8, label %_Z3fibi.exit.us, label %.lr.ph.i

_Z3fibi.exit.us:                                  ; preds = %_Z3fibi.exit.us, %4
  %i.01.us = phi i32 [ %11, %_Z3fibi.exit.us ], [ 0, %4 ]
  %10 = tail call i32 (i8*, ...)* @printf(i8* getelementptr inbounds ([15 x i8]* @.str1, i64 0, i64 0), i32 %7, double 1.000000e+00)
  %11 = add nsw i32 %i.01.us, 1
  %exitcond = icmp eq i32 %11, 10000
  br i1 %exitcond, label %.loopexit, label %_Z3fibi.exit.us

.lr.ph.i:                                         ; preds = %_Z3fibi.exit, %4
  %i.01 = phi i32 [ %16, %_Z3fibi.exit ], [ 0, %4 ]
  br label %12

; <label>:12                                      ; preds = %12, %.lr.ph.i
  %i.03.i = phi i32 [ 2, %.lr.ph.i ], [ %14, %12 ]
  %a.02.i = phi double [ 1.000000e+00, %.lr.ph.i ], [ %b.01.i, %12 ]
  %b.01.i = phi double [ 1.000000e+00, %.lr.ph.i ], [ %13, %12 ]
  %13 = fadd double %a.02.i, %b.01.i
  %14 = add nsw i32 %i.03.i, 1
  %exitcond.i = icmp eq i32 %14, %9
  br i1 %exitcond.i, label %_Z3fibi.exit, label %12

_Z3fibi.exit:                                     ; preds = %12
  %15 = tail call i32 (i8*, ...)* @printf(i8* getelementptr inbounds ([15 x i8]* @.str1, i64 0, i64 0), i32 %7, double %13)
  %16 = add nsw i32 %i.01, 1
  %exitcond2 = icmp eq i32 %16, 10000
  br i1 %exitcond2, label %.loopexit, label %.lr.ph.i

.loopexit:                                        ; preds = %_Z3fibi.exit, %_Z3fibi.exit.us, %2
  %.0 = phi i32 [ 1, %2 ], [ 0, %_Z3fibi.exit ], [ 0, %_Z3fibi.exit.us ]
  ret i32 %.0
}

declare i32 @puts(i8* nocapture) nounwind

declare i32 @atoi(i8* nocapture) nounwind readonly

declare i32 @printf(i8* nocapture, ...) nounwind

!0 = metadata !{metadata !"any pointer", metadata !1}
!1 = metadata !{metadata !"omnipotent char", metadata !2}
!2 = metadata !{metadata !"Simple C/C++ TBAA"}

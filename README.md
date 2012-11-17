auto-diagonalize
================

auto-diagonalize is an optimization pass written for LLVM.

Its purpose is to replace loops that can be modeled as linear dynamical
systems with a closed-form representation of the system.

The effect of this optimization is to transform simple _O(n)_ processes
into _O(lg(n))_ processes.

This is only possible if the _n x n_ transformation matrix of the system
has _n_ linearly independent eigenvectors.

### The Fibonacci Example

Take for example the iterative fibonacci algorithm;

```
; <label>:3                                       ; preds = %3, %.lr.ph
  %i.03 = phi i32 [ 2, %.lr.ph ], [ %5, %3 ]
  %a.02 = phi double [ 1.000000e+00, %.lr.ph ], [ %b.01, %3 ]
  %b.01 = phi double [ 1.000000e+00, %.lr.ph ], [ %4, %3 ]
  %4 = fadd double %a.02, %b.01
  %5 = add nsw i32 %i.03, 1
  %exitcond = icmp eq i32 %5, %2
  br i1 %exitcond, label %._crit_edge, label %3
```

The optimization pass eliminates this loop, replacing it with;

```
dgen:                                             ; preds = %.lr.ph
  %expt = uitofp i32 %n to double
  %eigvexpt = call double @llvm.pow.f64(double 0x3FF9E3779B97F4A6, double %expt)
  %pdn = fmul double 0x3FEB38880B4603E6, %eigvexpt
  %pdn1 = fmul double 0x3FE0D2CA0DA1530D, %eigvexpt
  %eigvexpt2 = call double @llvm.pow.f64(double 0xBFE3C6EF372FE94E, double %expt)
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
  br label %._crit_edge
```

Through constant folding and auto-vectorization, LLVM may produce something
even nicer.

### Limitations

This optimization currently only works on loops that are in [LCSSA]
(http://gcc.gnu.org/onlinedocs/gccint/LCSSA.html) form. It also only
handles loops which update floating point variables, though it should be
trivial to remove this restriction.

This pass does not attempt to handle transformation matrices which have
complex eigenvalues.

### Acknowledgements

This project makes heavy use of [Eigen]
(http://eigen.tuxfamily.org/index.php?title=Main_Page) and [LLVM](
http://llvm.org), which are both really fun to work with.

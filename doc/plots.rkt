#lang racket

(require plot)

(define xs (stream->list (in-range 1 200)))

(define (lg x)
  (/ (log x) (log 2)))

(define (mean pts)
  (/ (apply + pts) (length pts)))

(define (linear-regression pts)
  (define xbar (mean xs))
  (define ybar (mean pts))
  (define beta 
    (/ (foldl + 0 
              (map (lambda (i) (* (- i xbar)
                                  (- (list-ref pts (- i 1)) ybar)))
                   xs))
       (foldl + 0
              (map (lambda (i) (sqr (- i xbar))) xs))))
  (define alpha
    (- ybar (* beta xbar)))
  (printf "--- Linear Regression: y = alpha + beta*x --- ~n")
  (printf "alpha = ~a, beta = ~a ~n" alpha beta)
  (lambda (n)
    (+ alpha (* beta n))))

;; -O3, 50,000 repetitions per "n"
(define ys (list
            4102 1648 1083 1136 1491 1748 1380 1326 1387 1453 1496 1658 1690 1767 1757 1853 1994 1947 2068 2087 2152 2238 2296 2323 2400 2481 2516 2571 2652 2766 2777 2844 2901 3545 3174 3499 3689 3734 3798 3856 3912 3975 4051 4249 4197 4250 4301 4390 4437 4588 4582 4620 4681 4739 4823 4890 4921 4984 5120 5159 5294 5326 5336 5415 5460 5527 5601 5668 5711 5759 5811 5946 5951 6038 6104 6177 6237 6284 6369 6404 6482 6661 6604 6689 6714 6782 6867 6925 6995 7106 7094 7195 7232 7319 7371 7426 7470 7539 7708 7697 7730 7799 7878 7912 7978 8096 8241 8170 8232 8331 8355 8420 8487 8647 8632 8687 8808 8824 8865 8926 9087 9046 9143 9181 9308 9333 9361 9548 9518 9566 9660 9678 9562 9902 9900 9945 10017 10115 10119 10209 10421 10328 10358 10430 10556 10585 10715 10695 10818 10796 10919 10934 11118 11465 11156 11247 11277 11505 11394 11472 11506 11635 11650 11834 11762 11858 11875 11961 12111 12111 12160 12234 12261 12456 12419 12468 12540 12639 12747 12751 12788 12889 12933 12959 13168 13124 13245 13210 13446 13347 13432 13504 13575 15071 13644 14685 14630 13827 14081))

;; -fauto-diagonalize, 50,000 repetitions per "n"
(define zs (list
            3807 1029 1021 1077 1056 1048 1024 1021 1055 1043 1012 1492 1134 1047 1023 1041 1029 1067 1005 1027 1028 1013 1020 1035 1026 1012 1015 1049 1022 1015 1078 1053 1018 1007 1038 1037 1018 1023 1028 1069 1041 1033 1037 1023 1040 1013 1036 1023 1047 1017 1045 1015 1034 1007 1022 1019 1036 1041 1025 1064 1030 1014 1038 1025 1027 1051 1037 1048 1025 1034 1041 1061 1029 1003 1018 1023 1050 1051 1025 1024 1037 1028 1012 1016 1023 1038 1063 1035 1034 1018 1043 1047 1012 1044 1023 1032 1034 1077 1043 1030 1029 1031 1033 1036 1159 1030 1046 1008 1031 990 1036 1026 1030 1020 1002 1044 1021 1035 1022 1039 1036 1018 1053 1032 1052 995 1017 1022 1036 1026 1043 1027 1039 1073 1063 1061 1112 1044 1026 1064 1023 1031 1020 1014 1020 1035 1027 1053 1105 1033 1030 1023 1031 1017 1052 1024 1033 1047 1031 1007 1017 1032 1027 1043 1036 1007 1015 1069 1046 1017 1037 1037 1019 1021 1039 1020 1042 1024 1031 1018 1038 1138 1017 1039 1045 1050 1015 1037 1042 1047 1074 1024 1042 1032 1039 1041 1034 1038 1035))

;; -O3, int32_t fibi, 50,000 repetitions per "n"
(define as (list
            4400 1865 1065 1148 1179 1049 1187 1282 1275 1336 1376 1439 1449 1497 1466 1530 1634 1587 1610 1733 1738 1793 1789 1836 1909 1854 1915 1951 3103 4341 2083 2099 2197 2218 2244 2266 2334 2506 2433 2815 2882 2926 2946 2953 3076 3094 3075 3139 3194 3250 3280 2970 3358 3326 3497 3461 3521 3593 3618 3610 3641 3644 3792 3806 3830 3825 3886 3931 4005 4015 4096 4102 4083 4153 4164 4305 4361 4314 4425 4357 4413 4534 4566 4573 4618 4715 4701 4756 4774 5871 4860 4948 4920 4939 4973 5069 5028 5086 5248 5190 5222 5172 5412 5300 5439 5490 5390 5657 5546 5786 5682 5655 5711 5722 5725 5765 5860 5913 6054 5932 6223 5999 6099 6013 6100 6102 6145 6293 6194 6227 6292 6230 6289 6386 6546 6429 6452 6656 6594 6591 6652 6700 6659 6623 6642 6826 6785 6868 6833 6901 7024 7135 6921 6926 6968 7154 7344 7140 7167 7252 7302 7410 7235 7494 7379 7650 7616 7468 7485 7839 7609 7515 7756 7692 7972 7746 7790 7932 7800 8362 8247 8186 8066 7953 8282 8343 8137 8035 8475 8186 8172 8715 8587 8576 8526 8628 8547 8486 8694))

(define (plot-points pts title)
  (plot (list
         (points (map vector xs pts))
         (function (linear-regression pts) #:color "red"))
        #:y-min 0
        #:title title
        #:x-label "n"
        #:y-label "Microseconds Taken / 50,000 Calls"))

(plot-points ys "fib(n), clang++ -O3")
(plot-points zs "fib(n), clang++ -O3 -fauto-diagonalize")
(plot-points as "fibi(n), g++ -O3")
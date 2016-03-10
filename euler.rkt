#lang racket
(require redex)
(require redex/reduction-semantics)
(require "bmat-smat.rkt")
(require "junct-bmat.rkt")
(require (planet wmfarr/plt-linalg:1:13/matrix))
(provide (all-defined-out))

;integrates usng eulers method
(define (euler/ax+b a x b t dt ct)
  (cond
    ((>= ct t) x)
    (else (euler/ax+b a
                      (matrix-add x (matrix-scale (matrix-add (matrix-mul a x)b) dt))
                      b t dt (+ ct dt)))
    ))
;(mat->list (euler/ax+b (matrix 2 2 .5 0 0 1) (matrix 2 1 1 1) (matrix 2 1 0 0) 4 1 0))
;(mat->list (euler/ax+b (matrix 2 2 1 0 0 1) (matrix 2 1 1 1) (matrix 2 1 0 0) 4 .01 0))
;(exp 4)
;(mat->list (euler/ax+b (matrix 2 2 0 1 1 0) (matrix 2 1 0 0) (matrix 2 1 1 0) 4 1 0))

;integrates a junction using eulers
(define-metafunction bond
  euler/jct : junct any number number -> any
  ((euler/jct junct_start any_x number_t number_dt )
   ,(let ((ax+b (term (jct->ax+b junct_start))))
      (list (first ax+b) (euler/ax+b (second ax+b) (term any_x) (third ax+b) (term number_t) (term number_dt) 0)))
   ))
(define (disp/eul lst)
  (list (first lst) (mat->list (second lst))))
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (c 3)) ,(matrix 1 1 0) 1 1)))
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (c 3)) ,(matrix 1 1 0) 5 1)))
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (c 3)) ,(matrix 1 1 0) 10 1)))
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (c 3)) ,(matrix 1 1 0) 15 1)))
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (c 3)) ,(matrix 1 1 0) 20 1)))
;
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (c 3)) ,(matrix 1 1 0) 1 .01)))
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (c 3)) ,(matrix 1 1 0) 5 .01)))
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (c 3)) ,(matrix 1 1 0) 10 .01)))
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (c 3)) ,(matrix 1 1 0) 15 .01)))
;
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (0 (c 3))) ,(matrix 1 1 0) 20 1)));equivelent to above, works
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (0 (c 3) (r 5))) ,(matrix 1 1 0) 20 1)))
;(disp/eul (term (euler/jct (1 (v -7) (r 2) (1 (c 3))) ,(matrix 1 1 0) 20 1)));incorect result, can not have 1 conect to 1



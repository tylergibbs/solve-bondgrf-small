#lang racket
(require redex)
(require redex/reduction-semantics)
(require "junct-bmat.rkt")
(require (planet wmfarr/plt-linalg:1:13/matrix))
(provide (all-defined-out))


;makes a list of referencest to teh strings in bvec
(define (str#s bvec)
  (foldr (λ (x y z) (if (string? x) (cons y z) z))
         null
         bvec
         (build-list (length bvec) (λ (x) x))
         ))

;(str#s '("q" 0 0 -7 "p"))'(0 4)
;(str#s  '(-1 "q_2" 0 -7 0))'(2)

;makes a list of references to non0 numbers in bvec
(define (not0##s bvec)
  (foldr (λ (x y z) (if (not (or (string? x) (zero? x))) (cons (cons y x) z) z))
         null
         bvec
         (build-list (length bvec) (λ (x) x))
         ))
;(not0##s '("q" 0 0 -7 "p"))'(3)
;(not0##s  '(-1 "q_2" 0 -7 0));(0 3)

;gets the A matrix in Ax+B
(define (get-A bmat Arows smat sindex)
  (begin
    (map (λ (x y);x is ref to row in bmat, y is ref in smat
           (map (λ (z w);x is ref to col in bmat, y is ref in smat
                  (matrix-set! smat y w (matrix-ref bmat x z)))
                Arows
                sindex))
         Arows
         sindex;(1 2 3 ...)
         )
    smat))
;(mat->list (get-A (matrix 4 4
 ;                          1  2  3  4
  ;                         5  6  7  8
   ;                        9 10 11 12
    ;                      13 14 15 16) '(0 2) (matrix 2 2 0 0 0 0) '(0 1)))

;gets the B matrix in Ax+B
(define (get-B bmat Arows Brows/sclr smat sindex)
  (begin
    (map (λ (x);(cdr x) is the value of the source that created the row refed by (car x)
           (map (λ (y z);y is ref to col in bmat, z is ref in smat
                  (matrix-set! smat z 0 (+ (* (cdr x) (matrix-ref bmat y (car x))) (matrix-ref smat z 0))))
                Arows;only these colems added to Ax to get dx/dt since Ax is same rank
                sindex))
         Brows/sclr
         )
    smat))
;(mat->list (get-B (matrix 4 4
 ;                          1  2  3  4
  ;                         5  6  7  8
   ;                        9 10 11 12
    ;                      13 14 15 16) '(0 2) '((1 . 2)) (matrix 2 1 0 0)'(0 1)))

;gets a list of all strings in bvec
(define (get-svec bvec Arows)
  (map (λ (x) (list-ref bvec x))
       Arows
       ))
;(get-svec '("q" 0 0 -7 "p") '(0 4))'("q" "p")

;takes in the inverse of the big matrix of equations from bondgraph and their results and seperates
;them into A and B from Ax+B and names of each row
(define (bmat/bvec->ax+b bmat bvec)
  (let*((Arows (str#s bvec));each relevent row in bmat corisponds to a string in bvec(finding q' and p')
        (Brows/sclr (not0##s bvec));each value in B is from a source which is represented by a non 0 number in bvec
        (l (length Arows))
        (sindex (build-list (length Arows) (λ (x)x))));index of rows in resultent matrixes for convenience
  (list (get-svec bvec Arows)
        (get-A bmat Arows (make-matrix l l 0) sindex)
        (get-B bmat Arows Brows/sclr (make-matrix l 1 0) sindex))
    ))
(define (disp/ax+b lst)
  (list (first lst) (mat->list (second lst)) (mat->list (third lst))))
;
;(disp/ax+b (bmat/bvec->ax+b (matrix 4 4
;                           1  2  3  4
;                           5  6  7  8
;                           9 10 11 12
;                          13 14 15 16) '("q" -7 "p" 0)))
;;1VRC(v=-7 r=2 c=3)
;(mat->list (matrix-inverse (matrix 5 5
;                                   1 0 0 -1 0
;                                   0 0 0 0  2
;                                   0 0 1 1  0
;                                   0 3 0 1  0
;                                   0 0 0 1 -1)))
;(mat->list (matrix-inverse (matrix 4 4
;                                   0 0 0  2
;                                   0 1 1  0
;                                   3 0 1  0
;                                   0 0 1 -1)))
;(disp/ax+b (bmat/bvec->ax+b (matrix-inverse (matrix 5 5
;                                                  1 0 0 -1 0
;                                                  0 0 0 0  2
;                                                  0 0 1 1  0
;                                                  0 3 0 1  0
;                                                  0 0 0 1 -1))
;                                           '(0 "q" -7  0 0)
;                                           ))
;(disp/ax+b (bmat/bvec->ax+b (matrix-inverse (matrix 4 4 
;                                                   0 0 0  2
;                                                   0 1 1  0
;                                                   3 0 1  0
;                                                   0 0 1 -1))
;                                           '("q" -7  0 0)
;                                           ))

;combinds solve-bond and bmat/bvec->ax+b for convenience
(define-metafunction bond
  jct->ax+b : junct -> any
  ((jct->ax+b junct_start)
   ,(let ((bmat/bvec (term (solve-bond junct_start))))
      (bmat/bvec->ax+b
       (matrix-inverse (first bmat/bvec));must be inverted since in bmat/bvec is x'*bmat=x
       (vector->list (second bmat/bvec)))));changed to list for convenence
  )
;(disp/ax+b (term (jct->ax+b (1 (v -7) (c 3) (r 2)))))
;(disp/ax+b (term (jct->ax+b (1 (v -7) (l 3) (c 2) (r 5)))));same as book

#lang racket
(require redex)
(require redex/reduction-semantics)
(require (planet wmfarr/plt-linalg:1:13/matrix))
(provide (all-defined-out))

(define-language bond
  (junct (number 1port ... junct ...))
  (1port V F C L R)
  (V (v val))
  (F (f val))
  (C (c val))
  (L (l val))
  (R (r val))
  (val number)
  (pos number)
  )

;mat->list; matrix -> aloalon
(define (mat->list mat)
  (build-list (matrix-cols mat) (λ (m) (build-list (matrix-rows mat) (λ (n) (matrix-ref mat n m))))))

;if it is a 1 junct and contains a C port the C port must come first and
;the row of 1s will go at where the C port would have gone
;since C accumulates f into q and the matrix is ordered f e1 e2 ... for 1 junct
;same for 0 junct but with e and p
(define-metafunction bond
  same-CL : number 1port -> boolean
  ((same-CL 1 (c val_2)) #t)
  ((same-CL 0 (l val_2)) #t)
  ((same-CL number_1 1port_2) #f)
  )

;returns the position of the row of 1s for this junction(see same-CL)
(define-metafunction bond
  find-loc-1s : number (1port ...) -> pos
  ((find-loc-1s number_type (1port_1 ...))
   ,(let ((n (memf (λ (x) (term (same-CL  number_type ,x))) (term (1port_1 ...)))))
           ;n is location of first for whitch sameCL returns true
      (if n
          (- (+ 1 (length (term (1port_1 ...)))) (length n))
          0)))
  )
;(term (find-loc-1s 1 ((l 1) (r 2) (c 2))))
;(term (find-loc-1s 0 ((l 1) (r 2) (c 2))))
;(term (find-loc-1s 0 ((r 2) (c 2))))
;(term (find-loc-1s 1 ((r 2) (c 2))))

;encodes e0=e1+e2+e3/f0=f1+f2+f3 into the matrix
(define-metafunction bond
  add-row-1s : (pos pos ...) pos any -> any
  ((add-row-1s (pos_1port pos_rest ...) pos_111 any_mat)
  ,(void
      (map (λ (x) (matrix-set! (term any_mat) (term pos_111) x 1))
                (term (pos_rest ...))
                )
           (matrix-set! (term any_mat) (term pos_111) (term pos_1port) -1)
     )
  ))

;(define mat (matrix 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0))
;(term (add-row-1s (0 1 3) 2 ,mat))
;(mat->list mat)

;encodes equation for a 1port into the matrix
(define-metafunction bond
  add-1port-aux : 1port any any pos pos pos -> any
  ((add-1port-aux (v val_1) any_mat any_vec pos_lin pos_f pos_e);e1=V
   (,(matrix-set! (term any_mat) (term pos_lin) (term pos_e) 1)
    ,(vector-set! (term any_vec) (term pos_lin) (term val_1))
    )
   )
  ((add-1port-aux (f val_1) any_mat any_vec pos_lin pos_f pos_e);f1=F
   (,(matrix-set! (term any_mat) (term pos_lin) (term pos_f) 1)
    ,(vector-set! (term any_vec) (term pos_lin) (term val_1)))
   )
  ((add-1port-aux (c val_1) any_mat any_vec pos_lin pos_f pos_e);e*c=q_n
   (,(matrix-set! (term any_mat) (term pos_lin) (term pos_e) (term val_1))
   ,(vector-set! (term any_vec) (term pos_lin) (string-append "q_" (number->string (term pos_e)))))
   )
  ((add-1port-aux (l val_1) any_mat any_vec pos_lin pos_f pos_e);f*l=p_n
   (,(matrix-set! (term any_mat) (term pos_lin) (term pos_f) (term val_1))
    ,(vector-set! (term any_vec) (term pos_lin) (string-append "p_" (number->string (term pos_f)))))
   )
  ((add-1port-aux (r val_1) any_mat any_vec pos_lin pos_f pos_e);f*r-e=0
   (,(matrix-set! (term any_mat) (term pos_lin) (term pos_e) -1)
    ,(matrix-set! (term any_mat) (term pos_lin) (term pos_f) (term val_1))
    ,(vector-set! (term any_vec) (term pos_lin) 0)
    )
   )
  )
;(define mat1 (matrix 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ))
;(define vec1 (vector 0 0 0 0 0))
;(term (add-1port-aux (v -4) ,mat1 ,vec1 0 0 1))
;(mat->list mat1)
;(term (add-1port-aux (f -4) ,mat1 ,vec1 1 0 1))
;(mat->list mat1)
;(term (add-1port-aux (c 4) ,mat1 ,vec1 2 0 1))
;(mat->list mat1)
;(term (add-1port-aux (l 4) ,mat1 ,vec1 3 0 1))
;(mat->list mat1)
;(term (add-1port-aux (r 4) ,mat1 ,vec1 4 0 1))
;(mat->list mat1)
;vec1

;determins if pos_1 is a e or f position
;everything befor this step can be blind to type of junction
;however next step needs to know which location is the f and e
(define-metafunction bond
  add-1port : number 1port any any pos pos pos -> any
  ((add-1port 0 1port_1 any_mat any_vec pos_lin pos_1 pos_2)
   (add-1port-aux 1port_1 any_mat any_vec pos_lin pos_2 pos_1))
  ((add-1port 1 1port_1 any_mat any_vec pos_lin pos_1 pos_2)
   (add-1port-aux 1port_1 any_mat any_vec pos_lin pos_1 pos_2))
  )
;(define mat2 (matrix 2 2 0 0 0 0))
;(define vec2 (vector 0 0))
;(term (add-1port 0 (v -4) ,mat2 ,vec2 0 0 1))
;(mat->list mat2)
;vec2
;
;(define mat3 (matrix 2 2 0 0 0 0))
;(define vec3 (vector 0 0))
;(term (add-1port 1 (v -4) ,mat3 ,vec3 1 0 1))
;(mat->list mat3)
;vec3
;
;(define mat4 (matrix 2 2 0 0 0 0))
;(define vec4 (vector 0 0))
;(term (add-1port 0 (v -4) ,mat4 ,vec4 1 1 0))
;(mat->list mat4)
;vec4

;adds the 1port equations from a (single) junction into the matrix
(define-metafunction bond
  add-row-1ports : number (1port ...) any any pos pos pos -> any
  ((add-row-1ports number_type (1port_1 1port_2 ...) any_mat any_vec pos_111 pos_orglin pos_111)
   ;row 1s is in this pos, so add the 1port equation before the rest in reseved position
   ((add-1port number_type 1port_1 any_mat any_vec pos_orglin pos_orglin pos_111)
    (add-row-1ports number_type (1port_2 ...) any_mat any_vec pos_111 pos_orglin ,(add1 (term pos_111))))
         )
  ((add-row-1ports number_type (1port_1 1port_2 ...) any_mat any_vec pos_111 pos_orglin pos_curlin);add 1port eq in this poition
   ((add-1port number_type 1port_1 any_mat any_vec pos_curlin pos_orglin pos_curlin)
    (add-row-1ports number_type (1port_2 ...) any_mat any_vec pos_111 pos_orglin ,(add1 (term pos_curlin))))
   )
  ((add-row-1ports number_type () any_mat any_vec pos_111 pos_orglin pos_curlin);termination, no more 1ports to add
   ,(void))
  )
;(define mat5 (matrix 3 3 0 0 0 0 0 0 0 0 0))
;(define vec5 (vector 0 0 0))
;(term (add-row-1ports 1 ((v -4) (r 3)) ,mat5 ,vec5 0 0 1))
;(mat->list mat5)

;determins the location of the row of 1s, adds them, then executes adds the eq from the 1ports avoiding the row of 1s
(define-metafunction bond
  eval-1ports : (1port ...) pos number (pos ...) any any -> any
  ((eval-1ports (1port_1 ...) pos_lin number_type (pos_111 ...) any_mat any_vec)
   ,(let ((poft (term (find-loc-1s number_type (1port_1 ...)))));find location that row of ones should be reletive to begin row of jct
      (void (term (add-row-1s (pos_111 ...) ,(+ poft (term pos_lin)) any_mat));adds row of 1s to correct slot
             (term (add-row-1ports number_type (1port_1 ...) any_mat any_vec ,(+ poft (term pos_lin)) pos_lin ,(add1 (term pos_lin))))
             ))
   ))
;(define mat6 (matrix 3 3 0 0 0 0 0 0 0 0 0))
;(define vec6 (vector 0 0 0))
;(term (eval-1ports ((v -4)(r 3)) 0 1 (1 2) ,mat6 ,vec6))
;(mat->list mat6)

;recures solve over all the junction withing an origenal junct on the correct row number
(define-metafunction bond
  solve-aux : (junct ...) pos pos any any -> (pos ...)
  ((solve-aux (junct_1 junct_2 ...) pos_curlin pos_1add any_mat any_vec)
   ,(let ((nexlin (term (solve junct_1 pos_curlin pos_1add any_mat any_vec))))
      (cons (add1(term pos_curlin)) (term (solve-aux (junct_2 ...) ,nexlin pos_1add any_mat any_vec)))))
  ((solve-aux () pos_curlin pos_1add any_mat any_vec)
   ())
  )

;adds the equations from this jucntion and containd junctions to matrix,
;returns matrix line number to add next junction to(if applicable)
(define-metafunction bond
  solve : junct pos pos any any -> pos
  ((solve (number_type 1port_1 ... junct_1 ...) pos_lin pos_1add any_mat any_vec)
   ,(let* ((line# (add1 (term pos_lin)));starts adding equations at 1 line more that pos_lin
           (l (length (term (1port_1 ...))));# of 1ports in junct
           (pos-bonds-1port (build-list l (λ (x) (+ 1 line# x))));location of equations for 1 ports
           (pos-bonds (cons (term pos_1add);input(assumes that 1junct 0 junct alternate and this will be added not the same)
                            (append pos-bonds-1port;location of +ed variable for each 1port
                              (term (solve-aux (junct_1 ...) ,(+ line# l ) ,line# any_mat any_vec)))));location
                                ;of input for each junction off this junction, also adds there eqs into matirx
           ;used to place -1 1 1 1 .... row, since must add(have 1) each row exept first
           ;(it is the variable that is = for all 1ports in junct)
           (1port/void (term (eval-1ports (1port_1 ...) ,line#
                                    number_type ,pos-bonds any_mat any_vec)));adds in equations for 1ports and row of 1s
           )
      (last pos-bonds);last line edited in matrix by this solve
      )))
;(define mat7 (matrix 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0))
;(define vec7 (vector 0 0 0 0))
;(term (solve (1 (v -4) (r 3)) 0 0 ,mat7 ,vec7))
;(mat->list mat7)
;vec7
;
;(define mat8 (matrix 6 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0))
;(define vec8 (vector 0 0 0 0 0 0))
;(term (solve (1 (v -2) (c 4) (r 3) (l 5)) 0 0 ,mat8 ,vec8))
;(mat->list mat8)
;vec8
;
;(define mat9 (matrix 6 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0))
;(define vec9 (vector 0 0 0 0 0 0))
;(term (solve (0 (f -2) (r 3) (l 5) (c 4)) 0 0 ,mat9 ,vec9))
;(mat->list mat9)
;vec9

;used to find the size of matrix that equations are read into
;number of 1ports + number of junctions
(define-metafunction bond
  junct-size : junct -> number
  ((junct-size (number_1 1port_1 ... junct_1 ...))
   ,(foldr
     (λ (x y) (+ (term (junct-size ,x)) y))
     (+ 1 (length (term (1port_1 ...))))
     (term (junct_1 ...))))
  )
;(term (junct-size (1 (v 2) (r 4) (0 (l 3) (r 2)))))

;gets the corisponding big matrix from the bond graph contand in junct
(define-metafunction bond
  solve-bond : junct -> (any any)
  ((solve-bond junct_1)
   ,(let*((size (add1 (term (junct-size junct_1))));size of matrix is 1 larger becouse input to first junction is not stated(is 0)
          (mat (make-matrix size size 0))
          (vec (make-vector size 0))
          (s/v (term (solve junct_1 0 0 ,mat ,vec))))
      (begin (matrix-set! mat 0 0 1);makes matrix of full rank, is nullified since corisponds to 0 in vec when mat is inverted
   (list mat vec))))
  )
;conveniant display
(define (disp-jct x)
  (list (mat->list (first x)) (second x)))
;(disp-jct(term (solve-bond (1 (v -2) (c 4) (r 3) (l 5)))))
;(disp-jct(term (solve-bond (0 (v -2) (r 3) (l 5) (c 4)))))
;(disp-jct(term (solve-bond (1 (v -7) (r 2) (0 (r 3) (c 4))))))
;(disp-jct(term (solve-bond (0 (v -7) (1 (r 2) (c 3)) (1 (r 4) (l 5))))))
;(disp-jct(term (solve-bond (0 (v -7) (1 (r 2) (c 3)) (1 (r 4) (1 (c 5) (l 6)))))))

;(disp-jct(term (solve-bond (1 (v -7) (c 3) (r 2)))))
;(mat->list (matrix-inverse (first(term (solve-bond (1 (v -7) (l 3) (c 2) (r 5)))))))
;(disp-jct(term (solve-bond (1 (v -7) (l 3) (c 2) (r 5)))))
;~150 lines code
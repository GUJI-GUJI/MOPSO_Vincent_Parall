!------------------------------------------------------------------------------------------
! C This subroutine generate Gaussian noise series.
! C U: mean; G: standard variance; R: random seed, double precision; N: length of the series; 
! C A(:): output Gaussian noise series. 	
	SUBROUTINE NGRNS(U,G,R,N,A)
	DOUBLE PRECISION R,S,W,V,T
	INTEGER :: I,J,N,M
	REAL :: U,G
	REAL :: A(N)
	S=65536.0
	W=2053.0
	V=13849.0
	DO 20 J=1,N
	  T=0.0
	  DO 10 I=1,12
	    R=W*R+V
	    M=R/S
	    R=R-M*S
	    T=T+R/S
10	  CONTINUE
	  A(J)=U+G*(T-6.0)
20	CONTINUE
	
	RETURN
	END


!---------------------------------------------------------------
!This subroutine is to compute the inverse matrix of A
	SUBROUTINE BRINV(A,N,L)
	integer ::	I, J, K, L, N
	DIMENSION A(N,N)
	DOUBLE PRECISION A,T,D
	integer, allocatable :: IS(:), JS(:)
	allocate (IS(N))
	allocate (JS(N))
	L=1
	DO 100 K=1,N
	  D=0.0
	  DO 10 I=K,N
	  DO 10 J=K,N
	    IF (ABS(A(I,J)).GT.D) THEN
	      D=ABS(A(I,J))
	      IS(K)=I
	      JS(K)=J
	    END IF
10	  CONTINUE
	  IF (D+1.0.EQ.1.0) THEN
	    L=0
	    WRITE(*,20)
	    RETURN
	  END IF
20	  FORMAT(1X,'ERR**NOT INV')
	  DO 30 J=1,N
	    T=A(K,J)
	    A(K,J)=A(IS(K),J)
	    A(IS(K),J)=T
30	  CONTINUE
	  DO 40 I=1,N
	    T=A(I,K)
	    A(I,K)=A(I,JS(K))
	    A(I,JS(K))=T
40	  CONTINUE
	  A(K,K)=1/A(K,K)
	  DO 50 J=1,N
	    IF (J.NE.K) THEN
	      A(K,J)=A(K,J)*A(K,K)
	    END IF
50	  CONTINUE
	  DO 70 I=1,N
	    IF (I.NE.K) THEN
	      DO 60 J=1,N
	        IF (J.NE.K) THEN
	          A(I,J)=A(I,J)-A(I,K)*A(K,J)
	        END IF
60	      CONTINUE
	    END IF
70	  CONTINUE
	  DO 80 I=1,N
	    IF (I.NE.K) THEN
	      A(I,K)=-A(I,K)*A(K,K)
	    END IF
80	  CONTINUE
100	CONTINUE
	DO 130 K=N,1,-1
	  DO 110 J=1,N
	    T=A(K,J)
	    A(K,J)=A(JS(K),J)
	    A(JS(K),J)=T
110	  CONTINUE
	  DO 120 I=1,N
	    T=A(I,K)
	    A(I,K)=A(I,IS(K))
	    A(I,IS(K))=T
120	  CONTINUE
130	CONTINUE
	deallocate (IS)
	deallocate (JS)
	RETURN
	END   	

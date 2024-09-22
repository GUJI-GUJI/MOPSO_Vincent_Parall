MODULE ProblemBank
IMPLICIT NONE

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	This is the where MOPs are defined
SUBROUTINE MOP1_circles(fit, indv, viol, NF, NP)
	IMPLICIT NONE
	INTEGER, INTENT(IN)::NF, NP
	REAL, INTENT(IN):: indv(NP)
	REAL, INTENT(OUT):: fit(NF)
	INTEGER, INTENT(OUT)::viol
	
	fit(1) = indv(1)**2 + indv(2)**2
	fit(2) = (indv(1)-2.0)**2 + (indv(2)-2.0)**2
	viol = 0
END SUBROUTINE MOP1_circles

SUBROUTINE MOP2_Fonseca(fit, indv, viol, NF, NP)
	IMPLICIT NONE
	INTEGER, INTENT(IN)::NF, NP
	REAL, INTENT(IN):: indv(NP)
	REAL, INTENT(OUT):: fit(NF)
	INTEGER, INTENT(OUT)::viol
	
	fit(1) = 1-EXP(-((indv(1)- 1.0/(3.0)**0.5)**2.0 + &
					(indv(2) - 1.0/(3.0)**0.5)**2.0 + &
					(indv(3) - 1.0/(3.0)**0.5)**2.0));
	fit(2) = 1-EXP(-((indv(1)+ 1.0/(3.0)**0.5)**2.0 + &
					(indv(2) + 1.0/(3.0)**0.5)**2.0 + &
					(indv(3) + 1.0/(3.0)**0.5)**2.0));
	viol = 0
END SUBROUTINE MOP2_Fonseca

SUBROUTINE MOP3_Poloni(fit, indv, viol, NF, NP)
	IMPLICIT NONE
	INTEGER, INTENT(IN)::NF, NP
	REAL, INTENT(IN):: indv(NP)
	REAL, INTENT(OUT):: fit(NF)
	INTEGER, INTENT(OUT)::viol
	
	REAL:: 	A1 = 0.873648562314064, &	!= 0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2)
			A2 = 2.748572443268639, &	!=1.5*sin(1)-cos(1)+2*sin(2)-.5*cos(2)
			x, y ,B1, B2
	
	x = indv(1)
	y = indv(2)
	B1 = 0.5*SIN(x) - 2.0*COS(x) + SIN(y) - 1.5*COS(y)
	B2 = 1.5*SIN(x) - COS(x) + 2.0*SIN(y) - 0.5*COS(y)

	fit(1) = (1 + (A1-B1)**2 + (A2-B2)**2)
	fit(2) = ((x+3)**2 + (y+1)**2)
	
	viol = 0
END SUBROUTINE MOP3_Poloni

SUBROUTINE MOP4_Kursawe(fit, indv, viol, NF, NP)
	IMPLICIT NONE
	INTEGER, INTENT(IN)::NF, NP
	REAL, INTENT(IN):: indv(NP)
	REAL, INTENT(OUT):: fit(NF)
	INTEGER, INTENT(OUT)::viol
	
	REAL::a, b, x1, x2, x3
	
	a = 0.8
	b = 3.0
	x1 = indv(1)
	x2 = indv(2)
	x3 = indv(3)
	
	fit(1) = -10 * (exp(-0.2 * SQRT(x1**2.0 + x2**2.0)) &
    			  + exp(-0.2 * SQRT(x2**2.0 + x3**2.0)))
	fit(2) =(ABS(x1))**a + 5.0*SIN(x1**b) + &
    		(ABS(x2))**a + 5.0*SIN(x2**b) + &
    		(ABS(x3))**a + 5.0*SIN(x3**b)
    
	viol = 0
END SUBROUTINE MOP4_Kursawe

SUBROUTINE MOP5_Viennet(fit, indv, viol, NF, NP)
	IMPLICIT NONE
	INTEGER, INTENT(IN)::NF, NP
	REAL, INTENT(IN):: indv(NP)
	REAL, INTENT(OUT):: fit(NF)
	INTEGER, INTENT(OUT)::viol
	
	REAL:: x, y, temp

	y = indv(2)
	x = indv(1)
	temp = x**2+y**2
	
	fit(1) = 0.5 * temp + SIN(temp)
	fit(2) = (3.0*x-2.0*y+4.0)**2.0/8.0 + &
			 (x-y+1.0)**2.0/27.0 + 15.0
	fit(3) = 1.0 / (temp+1.0) - 1.1*EXP(-temp)

	viol = 0
END SUBROUTINE MOP5_Viennet


SUBROUTINE ZDT3_10D(fit, indv, viol, NF, NP)
	IMPLICIT NONE
	INTEGER, INTENT(IN)::NF, NP
	REAL, INTENT(IN):: indv(NP)
	REAL, INTENT(OUT):: fit(NF)
	INTEGER, INTENT(OUT)::viol
	
	REAL:: r
	INTEGER::dims = 10, i
	REAL, PARAMETER:: pi = 3.141592653589793


	fit(1) = indv(1)
	
	r = 0.0
	
	DO i = 2, dims
		r = r + indv(i)
	END DO
	
	fit(2) = (1.0+9.0*r/(dims-1)) * &
			(1-SQRT(indv(1)/(1.0+9.0*r/(dims-1))) - &
			(indv(1)/(1.0+9.0*r/(dims-1)))*SIN(10.0*pi*indv(1)));

	viol = 0
END SUBROUTINE ZDT3_10D

SUBROUTINE MOP_C4_Tanaka(fit, indv, viol, NF, NP)
	IMPLICIT NONE
	INTEGER, INTENT(IN)::NF, NP
	REAL, INTENT(IN):: indv(NP)
	REAL, INTENT(OUT):: fit(NF)
	INTEGER, INTENT(OUT)::viol
	
	REAL:: a, b
	
	a = 0.1
	b = 16.0

	fit(1) = indv(1)
	fit(2) = indv(2)

	viol = 0
	IF (- indv(1)**2.0 - indv(2)**2.0 + 1.0 + &
		a * COS(b*ATAN(indv(1)/indv(2)))>0.0) viol = 1
END SUBROUTINE 	MOP_C4_Tanaka

!SUBROUTINE MOP_18_Vincent(fit, indv, viol, NF, NP)
!	USE GLOBAL_VARIABLE
!    IMPLICIT NONE
!	INTEGER, INTENT(IN)::NF, NP
!	REAL, INTENT(IN):: indv(NP)
!	REAL, INTENT(OUT):: fit(NF)
!	INTEGER, INTENT(OUT)::viol
!	INTEGER::I
!    INTEGER:: indv_it(NP)
!    DO i=1,NP
!        Indv_it(i)=int(indv(i))   
!    END DO
!    Call MAIN_hydrodynamic(Indv_it,fit(1),fit(2))
!	viol = 0
!    DEALLOCATE(R_Upstream)
!    DEALLOCATE(R_Downstream)
!    DEALLOCATE(Upstream)
!	 DEALLOCATE(Downstream)
!    DEALLOCATE(Ele_Rel)
!    DEALLOCATE(Ca_cl)
!    DEALLOCATE(Arch_cl)
!    DEALLOCATE(Insi_cl)
!    DEALLOCATE(PumpST_cl)
!    DEALLOCATE(C_Opendegree)
!    DEALLOCATE(C_OutDischarge)
!    DEALLOCATE(Forebay_MAX,Forebay_MIN,Outlet_MAX,Outlet_MIN,DESIGN_ZIN)
!    DEALLOCATE(ST_STATE)
!END SUBROUTINE MOP_18_Vincent

!SUBROUTINE MOP_CG_Vincent(fit, indv, viol, NF, NP)
!	USE GLOBAL_VARIABLE
!    IMPLICIT NONE
!	INTEGER, INTENT(IN)::NF, NP
!	REAL, INTENT(IN):: indv(NP)
!	REAL, INTENT(OUT):: fit(NF)
!	INTEGER, INTENT(OUT)::viol
!	INTEGER::I
!    REAL:: indv_it(NP)
!    DO i=1,NP
!        Indv_it(i)=indv(i)   
!    END DO
!!    Call MAIN_hydrodynamic(Indv_it,fit(1),fit(2),fit(3))
!    !Call MAIN_hydrodynamic(Indv_it,fit(1),fit(2))
!    !这里的Indv_it有354个，即变量NP=354
!    !但进去就变成3个？
!	viol = 0
!!    IF(abs(Fit(1)-1000.0)<0.01)THEN
!!        viol=1    
!
!END SUBROUTINE MOP_CG_Vincent

SUBROUTINE MOP_TU(fit, indv, viol, NF, NP)
	USE GLOBAL_VARIABLE
    IMPLICIT NONE
	INTEGER, INTENT(IN)::NF, NP
	REAL, INTENT(IN):: indv(NP)
	REAL, INTENT(OUT):: fit(NF)
	INTEGER, INTENT(OUT)::viol
	INTEGER::I
    REAL:: indv_it(NP)
    DO i=1,NP
        Indv_it(i)=indv(i)
    END DO
	!Call MAIN_hydrodynamic(Indv_it,fit(1),fit(2),fit(3))
    Call Test_f(Indv_it,fit(1),fit(2),NP)
    !Call Test_f(Indv_it,fit(1))
    !这里的Indv_it有354个，即变量NP=354
    !但进去就变成3个？
	viol = 0
END SUBROUTINE MOP_TU

END MODULE ProblemBank





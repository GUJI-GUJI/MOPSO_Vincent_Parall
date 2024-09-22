!subroutine Rebuild_by_Time(Cal_Time)
!    Boundary_in.End_Time=Boundary_in.Start_Time+Cal_Time
!    
!    OPEN(UNIT=11,FILE="INPUTFILE/BOUNDARY.txt")
!        READ(11,*),READ(11,*),READ(11,*),READ(11,*),READ(11,*),READ(11,*),READ(11,*),READ(11,*)
!        CALL Boundary_Auto_Update(Boundary_in,"NUM")
!        DO i=1,Boundary_in.NUM !¶Á±ß½çÌõ¼þ
!            READ(11,*)Up_tmp,Down_tmp
!            IF(Boundary_in.Upsign==0)THEN
!                Boundary_in.UpQ(i)=Up_tmp
!            ELSE IF(Boundary_in.Upsign==1)THEN
!                Boundary_in.UpZ(i)=Up_tmp
!            END IF
!            IF(Boundary_in.Downsign==0.or.Boundary_in.Downsign==2)THEN
!                Boundary_in.DownQ(i)=Down_tmp
!            ELSE IF(Boundary_in.Downsign==1)THEN
!                Boundary_in.DownZ(i)=Down_tmp
!            END IF
!        END DO
!    CLOSE(11)
!    CALL Model_Meshing
!    
!    
!    
!    
!    
!end subroutine
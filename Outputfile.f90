    !SUBROUTINE Outputfile(Initial_Result_main,Com_element_main,Dynamic_Result_out)
    !USE GLOBAL_VARIABLE
    !
    !IMPLICIT NONE
    !TYPE(Initial_Result)::Initial_Result_main
    !TYPE(Dynamic_Result)::Dynamic_Result_out
    !TYPE(Com_element)::Com_element_main
    !INTEGER::i,k,mm
    !CHARACTER(LEN=100),allocatable :: Q_char(:),Z_char(:)
    !CHARACTER(LEN=100):: temp
    !
    !!!******************输出初始计算水动力过程********************
    !!OPEN(UNIT=20,FILE="OUTPUTFILE/InitialRESULT.txt")
    !!    WRITE(20,"(A10,3X,A10,3X,A7,3X,A7,3X,A7,3X,A7,3X,A7,3X,A7,3X)")"断面序号","里程","水深","流量","流速","水位","总水头","蓄量"
    !!    DO i=1,Com_element_main.Element_NUM+1
    !!        WRITE(20,"(I6,7X,F10.3,3X,F7.3,3X,F7.3,3X,F7.3,3X,F7.3,3X,F7.3,3X,F15.3,3X)")i,&
    !!            Initial_Result_main.L(i),&
    !!            Initial_Result_main.h(i),&
    !!            Initial_Result_main.Q(i),&
    !!            Initial_Result_main.vc(i),&
    !!            Initial_Result_main.Z(i),&
    !!            Initial_Result_main.Z(i)+Initial_Result_main.vc(i)**2/2.0/9.81,&
    !!            Initial_Result_main.Vol(i)
    !!    END DO
    !!CLOSE(20)
    !
    !!******************输出初始计算水动力过程********************
    !OPEN(UNIT=20,FILE="OUTPUTFILE/InitialRESULT.csv")
    !WRITE(20,'(*(G0.5,:,","))')"断面序号","里程","水深","流量","流速","水位","总水头","蓄量"
    !DO i=1,Com_element_main.Element_NUM+1
    !    WRITE(20,'(*(G0.5,:,","))')i,&
    !        Initial_Result_main.L(i),&
    !        Initial_Result_main.h(i),&
    !        Initial_Result_main.Q(i),&
    !        Initial_Result_main.vc(i),&
    !        Initial_Result_main.Z(i),&
    !        Initial_Result_main.Z(i)+Initial_Result_main.vc(i)**2/2.0/9.81,&
    !        Initial_Result_main.Vol(i)
    !END DO
    !CLOSE(20)
    !
    !
    !OPEN(UNIT=44,FILE="OUTPUTFILE/InitialVolume.txt")
    !DO i=1,Com_element_main.Controlgate_NUM+1
    !    WRITE(44,"(I10,3X,F12.2,3X)")i,Initial_Result_main.Pool_Volume(i)
    !END DO
    !CLOSE(44)
    !
    !!!******************输出非恒定流计算水动力过程******************
    !!OPEN(UNIT=28,FILE="OUTPUTFILE/DynamicRESULT.txt")
    !!DO k=1,Dynamic_Result_out.NUM,1
    !!    IF(MOD((k-1)*Dynamic_Result_out.DeltaT,Out_DeltaT)==0)THEN
    !!        IF(Output_sign==0)THEN
    !!            WRITE(28,"(f7.2,7x,60(f7.3,3x,f7.3,3x,f7.3,7x))")(k-1)*Dynamic_Result_out.DeltaT/REAL(Out_DeltaT),(Dynamic_Result_out.Zup(i,k),Dynamic_Result_out.Zdown(i,k),Dynamic_Result_out.Q(i,k),i=1,Com_element_main.Controlgate_NUM),Dynamic_Result_out.Zbd(1,k),Dynamic_Result_out.Qbd(1,k)
    !!        ELSEIF(Output_sign==1)THEN
    !!            WRITE(28,"(f7.2,7x,1000(f7.3,3x,f7.3,7x))")(k-1)*Dynamic_Result_out.DeltaT/REAL(Out_DeltaT),(Dynamic_Result_out.ZZ(i,k),Dynamic_Result_out.QQ(i,k),i=1,Com_element_main.Element_NUM-Com_element_main.Bridge_NUM),Dynamic_Result_out.Zbd(1,k),Dynamic_Result_out.Qbd(1,k)
    !!        END IF
    !!    END IF
    !!END DO
    !!CLOSE(28)
    !
    !!******************输出非恒定流计算水动力过程******************
    !OPEN(UNIT=28,FILE="OUTPUTFILE/DynamicRESULT.csv")
    !ALLOCATE (Q_char(Com_element_main.Element_NUM-Com_element_main.Bridge_NUM+1))
    !ALLOCATE (Z_char(Com_element_main.Element_NUM-Com_element_main.Bridge_NUM+1))
    !DO i = 1,Com_element_main.Element_NUM-Com_element_main.Bridge_NUM+1!
    !    write(temp,'(i0)') i
    !    Q_char(i)='Q'//temp
    !    Z_char(i)='Z'//temp
    !END DO
    !if(Output_sign==1)then
    !    WRITE(28,'(*(G0.5,:,","))')'时间序列',(Z_char(i),Q_char(i),i=1,Com_element_main.Element_NUM-Com_element_main.Bridge_NUM+1)
    !end if
    !
    !DO k=1,Dynamic_Result_out.NUM,1
    !    IF(MOD((k-1)*Dynamic_Result_out.DeltaT,Out_DeltaT)==0)THEN
    !        IF(Output_sign==0)THEN
    !            WRITE(28,'(*(G0.5,:,","))')(k-1)*Dynamic_Result_out.DeltaT/REAL(Out_DeltaT),(Dynamic_Result_out.Zup(i,k),Dynamic_Result_out.Zdown(i,k),Dynamic_Result_out.Q(i,k),i=1,Com_element_main.Controlgate_NUM),Dynamic_Result_out.Zbd(1,k),Dynamic_Result_out.Qbd(1,k)
    !        ELSEIF(Output_sign==1)THEN
    !            WRITE(28,'(*(G0.5,:,","))')(k-1)*Dynamic_Result_out.DeltaT/REAL(Out_DeltaT),(Dynamic_Result_out.ZZ(i,k),Dynamic_Result_out.QQ(i,k),i=1,Com_element_main.Element_NUM-Com_element_main.Bridge_NUM),Dynamic_Result_out.Zbd(1,k),Dynamic_Result_out.Qbd(1,k)
    !        END IF
    !    END IF
    !END DO
    !CLOSE(28)
    !
    !
    !
    !!******************输出各渠池蓄量过程******************
    !OPEN(UNIT=31,FILE="OUTPUTFILE/VolumeRESULT.txt")
    !DO k=1,Dynamic_Result_out.NUM,1
    !    IF(MOD((k-1)*Dynamic_Result_out.DeltaT,Out_DeltaT)==0)THEN
    !        WRITE(31,"(f7.2,7x,70(f15.2,3x))")(k-1)*Dynamic_Result_out.DeltaT/REAL(Out_DeltaT),(Dynamic_Result_out.Pool_Volume(i,k),i=1,Com_element_main.Controlgate_NUM+1),Dynamic_Result_out.Total_Volume(k)
    !    END IF
    !END DO
    !CLOSE(31)
    !
    !
    !END

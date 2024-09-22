    !Subroutine MAIN_hydrodynamic(Bav,D_Target,T_Target,M_Target)
    Subroutine MAIN_hydrodynamic(Bav,D_Target,M_Target)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    !/////һά����ˮ���������ˮ����ģ�⼰�����ʶ�///////
    INTEGER::i,j,k,gg,fpDynamicResult,fpControlGate,fpQControlGate,fpVolume,&
        fpBestDynamicResult,fpBestControlGate,fpBestQControlGate,fpBestVolume
    REAL::A_n,B_n
    REAL::D_Target  !Ŀ�꺯��1����Ŀ��ˮλ�����ƫ����С
    REAL::M_Target  !Ŀ�꺯��2����Ŀ��ˮλ��ƽ��ƫ����С
    !REAL::T_Target   !Ŀ�꺯��3��բ�ŵ��ش�������
    REAL::Bav(1:Controlgate_NUM*Controltime)!Controltime���ʱ�䲽�� 3secs
    !��ˮ���������ȫ�ֱ�����ת��Ϊ�ֲ�����
    REAL,ALLOCATABLE::GAV_Opendegree(:) !�����������
    REAL,ALLOCATABLE::Deltadegree(:,:)
    REAL,ALLOCATABLE::Opendegree(:,:)   !բ�ſ��ȹ���
    REAL,ALLOCATABLE::C_Opendegree(:,:) !�������ʹ��բ�ſ���
    REAL,ALLOCATABLE::ST_STATE(:,:)     !�洢������բÿСʱբǰˮλ��բ��ˮλ
    REAL,ALLOCATABLE::DS_STATE(:)       !�洢�����α߽�ˮλ
    REAL,ALLOCATABLE::Pool_Volume(:,:)
    REAL::QControlGate((End_Time-Start_Time)/3600+1,2*Controlgate_NUM) !բ����������
    !��¼����ˮ�����ɹ���
    fpDynamicResult=50
    fpBestDynamicResult=150
    fpControlGate=51
    fpBestControlGate=151
    fpQControlGate=52
    fpBestQControlGate=152
    fpVolume=53
    fpBestVolume=153
    OPEN(UNIT = fpDynamicResult, FILE = ('out/Dynamicresult.txt'))
    OPEN(UNIT = fpControlGate, FILE = ('out/ControlGate.txt'))     !������բբǰբ��ˮλ
    OPEN(UNIT = fpQControlGate, FILE = ('out/QControlGate.txt'))   !������բ����
    OPEN(UNIT = fpVolume, FILE = ('out/Volume.txt'))               !����������
    OPEN(UNIT = fpBestDynamicResult, FILE = ('out/BestDynamicresult.txt'))
    OPEN(UNIT = fpBestControlGate, FILE = ('out/BestControlGate.txt'))     !��Ѹ�����բբǰբ��ˮλ
    OPEN(UNIT = fpBestQControlGate, FILE = ('out/BestQControlGate.txt'))   !��Ѹ�����բ����
    OPEN(UNIT = fpBestVolume, FILE = ('out/BestVolume.txt'))               !����������
    !*****************************************��
    !*****�����ع��̴��ݵ�����ˮ��������******��
    !*****************************************��
    ALLOCATE(GAV_Opendegree(Controlgate_NUM*Controltime))
    !    write(*,*)D_Target,M_Target
    DO i=1,Controlgate_NUM*Controltime
        GAV_Opendegree(i)=Bav(i)
    END DO
    ALLOCATE(Deltadegree(Controlgate_NUM,Boundary_NUM))
    ALLOCATE(Opendegree(Controlgate_NUM,Boundary_NUM))
    DeltaDegree(:,:)=0.0  !��ʼ��ֵΪ0
    
    !��һ����ʲô��˼��
    !���������ļ�����բֻ��Ҫ��һ�����ȣ�����Ļ��Զ����룿
    DO i=1,Controlgate_NUM
        DO j=1,Controltime    !����6�ε�������
            DeltaDegree(i,4*j-3)=GAV_Opendegree(Controltime*(i-1)+j)
            DeltaDegree(i,4*j-2)=0.0
            DeltaDegree(i,4*j-1)=0.0
            DeltaDegree(i,4*j-0)=0.0
        END DO
    END DO

    Opendegree(:,1)=Opendegree_INIT(:,1)!����ĳ�ʼ����
    DO i=1,Controlgate_NUM
        Opendegree(i,1)=Opendegree(i,1)+DeltaDegree(i,1)
        !���������м�������բ�����������м�������ʱ��
        DO j=2,Boundary_NUM
            Opendegree(i,j)=Opendegree(i,j-1)+DeltaDegree(i,j)
        ENDDO
    ENDDO
    !////����۲�һ��
    ! WRITE(*,*)"��ǰ�������ӿ��ȵ��ع��̣�����"
    !DO j=1,Boundary_NUM
    !WRITE(*,"(13(F7.3,2X))")(Opendegree(i,j),i=1,Controlgate_NUM)
    !ENDDO

    !//////բ�ع��̼��ܣ�Ӧ�÷���ˮ����������
    IF(Boundary_Step/=DeltaT)THEN
        ALLOCATE(C_Opendegree(Controlgate_NUM,Calculate_NUM))
        k=Boundary_Step/DeltaT
        DO j=1,Calculate_NUM-2
            C_Opendegree(:,j)=Opendegree(:,j/k+1)+(MOD(j,k)-1)*(Opendegree(:,j/k+2)-Opendegree(:,j/k+1))/k
        END DO
        C_Opendegree(:,Calculate_NUM-1)=Opendegree(:,Boundary_NUM)-(Opendegree(:,Boundary_NUM)-Opendegree(:,Boundary_NUM-1))/k
        C_Opendegree(:,Calculate_NUM)=Opendegree(:,Boundary_NUM)
    ELSEIF(Boundary_Step==DeltaT)THEN
        ALLOCATE(C_Opendegree(Controlgate_NUM,Calculate_NUM))
        k=Boundary_Step/DeltaT
        DO j=1,Calculate_NUM
            C_Opendegree(:,j)=Opendegree(:,j)
        END DO
    END IF
    !/////////////////////////////////////////////////////////////////////////!
    !/////////����۲�һ��//////////
    !   WRITE(*,*)"��ʼ�Ǻ㶨�����㣡����"
    !******************************!
    !*********�Ǻ㶨������*********��
    !******************************!
    
    ALLOCATE(ST_STATE((End_Time-Start_Time)/3600+1,2*Controlgate_NUM))   !�洢ÿСʱ������բ��բǰˮλ��բ��ˮλ�����ڼ���Ŀ�꺯��
    ALLOCATE(DS_STATE((End_Time-Start_Time)/3600+1))    !�洢���α߽�ˮλ�����ڼ���Ŀ�꺯��
    
    ALLOCATE(Pool_Volume(Controlgate_NUM+1,Calculate_NUM))
    !ALLOCATE(QControlGate((End_Time-Start_Time)/3600+1,2*Controlgate_NUM))
    CALL Dynamic_Cal(C_Opendegree,ST_STATE,DS_STATE,Pool_Volume,QControlGate)!���˿��ȣ�����ȫ�����ֵ
    !/////�Ǻ㶨������


    !******************************!
    !*********Ŀ�꺯������**********��
    !******************************!
    CALL Function_Target(GAV_Opendegree,ST_STATE,DS_STATE,D_Target,M_Target)!,QControlGate



    IF(D_Target<=MIN_Target1)THEN
        write(*,*)D_Target,M_Target
        BestST_STATE=ST_STATE
        BestOpendegree=Opendegree
        BestQControlGate=QControlGate
        BestVolume=Pool_Volume
        write(*,*)"�Ѹ�������ֵ"
    END IF

    !*********����Ǻ㶨�����********��
    If(endflag==1)then
        !////���ˮλ�仯////
        write(fpDynamicResult,*)"----"
        write(fpDynamicResult,*)D_Target,M_Target
        CALL writeTOfile(fpDynamicResult,ST_STATE,(End_Time-Start_Time)/3600+1,2*Controlgate_NUM)
        !////���բ�ű仯////
        write(fpControlGate,*)"----"
        write(fpControlGate,*)D_Target,M_Target
        DO j=1,Boundary_NUM
            WRITE(fpControlGate,"(59(F7.3,2X))")(Opendegree(i,j),i=1,Controlgate_NUM)
        ENDDO
        !////���բ�������仯////
        write(fpQControlGate,*)"----"
        write(fpQControlGate,*)D_Target,M_Target
        CALL writeTOfile(fpQControlGate,QControlGate,(End_Time-Start_Time)/3600+1,2*Controlgate_NUM)
        !////������������仯////
        write(fpVolume,*)"----"
        write(fpVolume,*)D_Target,M_Target
        CALL writeTOfile(fpVolume,Pool_Volume,Controlgate_NUM+1,Calculate_NUM)

        IF(bestoutputflag==1)THEN
            bestoutputflag=0
            !////������ˮλ�仯////
            CALL writeTOfile(fpBestDynamicResult,BestST_STATE,(End_Time-Start_Time)/3600+1,2*Controlgate_NUM)
            !////������բ�ű仯////
            DO j=1,Boundary_NUM
                WRITE(fpBestControlGate,"(59(F7.3,2X))")(BestOpendegree(i,j),i=1,Controlgate_NUM)
            ENDDO
            !////������բ�������仯////
            CALL writeTOfile(fpBestQControlGate,BestQControlGate,(End_Time-Start_Time)/3600+1,2*Controlgate_NUM)
            !////������բ�������仯////
            CALL writeTOfile(fpBestVolume,BestVolume,Controlgate_NUM+1,Calculate_NUM)
        ENDIF
    end if


    DEALLOCATE(C_Opendegree)
    DEALLOCATE(ST_STATE)
    DEALLOCATE(DS_STATE)
    DEALLOCATE(GAV_Opendegree)
    DEALLOCATE(Opendegree)
    DEALLOCATE(Deltadegree)
    DEALLOCATE(Dynamic_Result_Zup)
    DEALLOCATE(Dynamic_Result_Zdown)
    DEALLOCATE(Dynamic_Result_Q)
    DEALLOCATE(Dynamic_Result_Zbd)
    DEALLOCATE(Dynamic_Result_Qbd)
    DEALLOCATE(Dynamic_Result_ZZ)
    DEALLOCATE(Dynamic_Result_QQ)
    END
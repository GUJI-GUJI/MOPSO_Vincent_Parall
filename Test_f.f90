    Subroutine Test_f(Input,out_lenth,out_time,NP)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    INTEGER::i,j,k
    REAL::Input(NP)
    real,allocatable::OldOpendegree(:,:)
    INTEGER::NP

    real::pollu_Lenth,arrive_time!��Ⱦ������
    REAL::out_lenth,out_time

    !REAL::T1_Gain  !Ŀ�꺯��1���������
    !REAL::x(3)     !x1,x2,x3�����
    !real::y1
    TYPE(Boundary)::Boundary_temp

    INTERFACE
    !�Բ�ͬ��С��һά����������Բ�ֵ
    SUBROUTINE Linear_Interpolation(array1,array2)
    REAL,ALLOCATABLE::array1(:),array2(:)
    END SUBROUTINE
    !�Բ�ͬ��С�Ķ�ά����������Բ�ֵ
    SUBROUTINE Linear_Interpolation_2d(array1,array2,direction)
    REAL,ALLOCATABLE::array1(:,:),array2(:,:)
    INTEGER::direction
    END SUBROUTINE
    END INTERFACE
    
    ALLOCATE(OldOpendegree(Controlgate_NUMber,Boudarys_Number))!δ��ֵ����
    ALLOCATE(Boundary_temp.Opendegree(Controlgate_NUMber,Boundary_cal.NUM))


    OldOpendegree(:,1)=Opendegree_INIT(:)!����ĳ�ʼ����

    Do i=1,Controlgate_NUMber
        OldOpendegree(i,1)=OldOpendegree(i,1)+input((i-1)*Boudarys_Number+1)
    ENDDO
    !���������м�������բ�����������м�������ʱ��
    Do i=1,Controlgate_NUMber
        DO j=2,NP
            OldOpendegree(i,j)=OldOpendegree(i,j-1)+input((i-1)*Boudarys_Number+j)
            IF(OldOpendegree(i,j)<0) THEN
                OldOpendegree(i,j)=0
            ENDIF
        END DO
    ENDDO
     OldOpendegree(:,NP+1:)=OldOpendegree(1,NP)
    
    !��ֵ
    k=Boundary_in.step/Boundary_cal.step
    do i = 1,Boundary_cal.NUM-2
        Boundary_temp.Opendegree(:,i)=OldOpendegree(:,i/k+1)+(MOD(i,k)-1)*(OldOpendegree(:,i/k+2)-OldOpendegree(:,i/k+1))/k
    end do
    Boundary_temp.Opendegree(:,Boundary_cal.NUM-1)=OldOpendegree(:,Boudarys_Number)-(OldOpendegree(:,Boudarys_Number)-OldOpendegree(:,Boudarys_Number-1))/k
    Boundary_temp.Opendegree(:,Boundary_cal.NUM)=OldOpendegree(:,Boudarys_Number)


    !write(41,'(*(G0.5,:,","))') (1000*Boundary_temp.Opendegree(1,j),j=1,Boundary_cal.NUM,6)!���ÿһ��բ�ſ���

    pollu_Lenth=0    !�洢���α߽�ˮλ�����ڼ���Ŀ�꺯��
    CALL Dynamic_Cal(Boundary_temp,pollu_Lenth,arrive_time)!���˿��ȣ�����ȫ�����ֵ
    out_lenth=pollu_Lenth
    out_time=arrive_time
    IF(pollu_Lenth<10000) startflag=0
    IF(pollu_Lenth<MIN_Target1) MIN_Target1=pollu_Lenth
    IF(arrive_time<MIN_Target2) MIN_Target2=arrive_time

    !y1=(x(1)-1)**2+(x(2)-2)**2+(x(3)-3)**2+14
    !IF(y1<100) startflag=0
    !IF(y1<MIN_Target1) MIN_Target1=y1
    !T1_Gain=y1

    END Subroutine


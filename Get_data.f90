    subroutine Get_data(Length_Time,Out_dOpendegree,NF,NP,length,outOpendegree)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    INTEGER, INTENT(IN)::NF,NP,length
    REAL, INTENT(IN)::Length_Time(NF*length)
    REAL, INTENT(IN)::Out_dOpendegree(NP*length)

    real,allocatable::polluLength(:)
    real,allocatable::time(:)
    real,allocatable::dOpendegree(:,:)
    real,allocatable::Opendegree(:,:)
    real::outOpendegree(NP)

    integer,allocatable::sortIndex(:)

    INTEGER:: i,j,n,k
    
    real::temp
    integer::temp2
    real,allocatable::temp3(:)



    Allocate(polluLength(length))
    Allocate(time(length))
    Allocate(dOpendegree(length,NP))
    Allocate(Opendegree(length,NP))


    do i = 1,length
        polluLength(i)=Length_Time(2*i-1)
        time(i)=Length_Time(2*i)
    end do

    do i = 1,length
        dOpendegree(i,:)=Out_dOpendegree((i-1)*NP+1:i*NP)
    end do

    Allocate(sortIndex(length))
    Allocate(temp3(NP))

    do i=1,length
        sortIndex(i)=i
    end do

    do i=1,length-1
        do j=i+1,length
            if (polluLength(i).gt.polluLength(j)) then
                temp = polluLength(i)
                polluLength(i) = polluLength(j)
                polluLength(j) = temp

                temp=time(i)
                time(i) = time(j)
                time(j) = temp

                temp3(:) = dOpendegree(i,:)
                dOpendegree(i,:) = dOpendegree(j,:)
                dOpendegree(j,:) = temp3(:)

                temp2=sortIndex(i)
                sortIndex(i) = sortIndex(j)
                sortIndex(j) = temp2
            endif
        enddo
    enddo


    DO i=1,Com_element_main.Controlgate_NUM
        Opendegree(:,1+(i-1)*(NP/Com_element_main.Controlgate_NUM))=Opendegree_INIT(i)+dOpendegree(:,1+(i-1)*(NP/Com_element_main.Controlgate_NUM))
    ENDDO


    do j=2,NP/Com_element_main.Controlgate_NUM
        DO i=1,Com_element_main.Controlgate_NUM
            Opendegree(:,j+(i-1)*(NP/Com_element_main.Controlgate_NUM))=Opendegree(:,j-1+(i-1)*(NP/Com_element_main.Controlgate_NUM))+dOpendegree(:,j+(i-1)*(NP/Com_element_main.Controlgate_NUM))
        ENDDO
    end do


    DO i=1,size(Opendegree,1)
        DO j=1,size(Opendegree,2)
            IF(Opendegree(i,j)<0)THEN
                Opendegree(i,j)=0
            ENDIF
        ENDDO
    ENDDO

    !ÿ25mins����һ��
    !����������ʱ�����Ⱦ��ɢȨ�أ�����1����բ�Ź�С����ȫ������Ⱦ��ɢ��Χ��С
    !����100�������ȫ���ǵ���ʱ��
    DO i=1,size(Opendegree,1)
        do j=1,5
            DO k=1,Com_element_main.Controlgate_NUM
                Opendegree(i,(j-1)*5+2+(k-1)*NP/Com_element_main.Controlgate_NUM  :(j-1)*5+5+(k-1)*NP/Com_element_main.Controlgate_NUM)=Opendegree(i,(j-1)*5+1+(k-1)*NP/Com_element_main.Controlgate_NUM)
            ENDDO
        end do
    ENDDO
    
    !open(unit=71,file="outputfile/pareto-frontier.csv")
    !write(71,*)"�������",",","��Ⱦ����ɢ����(m)",",","��Ⱦ�ŵ���ʱ��(mins)",",","�������С�����������>"
    !do i=1,length
    !    WRITE(71,'(*(G0.5,:,","))') i,polluLength(i),Time(i),(Opendegree(i,j),j=1,NP)
    !end do
    !close(71)
    
    open(unit=71,file="outputfile/selected-condition.txt")
    write(71,*)"��Ⱦ����ɢ��Χ��m��"
    write(71,*)polluLength(arrivetime_weight)
    write(71,*)"��Ⱦ�ﵽ��ʱ�䣨min��"
    write(71,*)Time(arrivetime_weight)
    close(71)
    
    outOpendegree=Opendegree(arrivetime_weight,:)
    
    !open(unit=71,file="outputfile/GateOpenDegree.csv")
    !write(71,*)"����բ����"
    !do i=1,length
    !    WRITE(71,'(*(G0.5,:,","))') i,outOpendegree(i)
    !end do
    !close(71)
    
    
    
    
    
    
    end subroutine
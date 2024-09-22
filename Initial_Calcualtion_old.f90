SUBROUTINE Initial_Cal  !wjb0413
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    
    INTEGER::i,j,k,gg,mm
    !/////��ʼ����
    !��������ˮλ����������洢�ռ�
    ALLOCATE(h(Element_NUM+1))
    ALLOCATE(Q(Element_NUM+1))
    ALLOCATE(Z(Element_NUM+1))
    ALLOCATE(vc(Element_NUM+1))
    ALLOCATE(L(Element_NUM+1))
    ALLOCATE(Vol(Element_NUM+1))
    !����Ԫ���洢��������洢�ռ�
    ALLOCATE(nnn(Element_NUM))
	Waterdepth_flag=0
    !/////��ʼ��������������
    IF(Boundary_sign_UP==0)THEN    !�������α߽�Ĳ�ͬѡ��������ֵ
        Q(1)=Upstream(1)
    ELSE
        Q(1)=Initial_Discharge
    END IF
    j=0
    DO i=2,Element_NUM
        IF(Ele_Rel(i,2)==4)THEN
            j=j+1
            Q(i)=Q(i-1)-OutDischarge(j,1)
        Else
            Q(i)=Q(i-1)
        END IF
    END DO
    Q(Element_NUM+1)=Q(Element_NUM)

    !//////��ʼˮλ�ļ��㣬�Լ������εķֶΣ�Ŀ����ȷ���������ˮλ
    !//////ˮλ�ļ����Ǵ���������������
    !//////����ˮλ�߽�ֱ�Ӹ�ֵˮλ
    h(Element_NUM+1)=D_elevation
    Waterdepth_flag=Waterdepth_flag+1
    Waterdepth(Waterdepth_flag,1)=h(Element_NUM+1)
    Waterdepth(Waterdepth_flag,2)=Q(Element_NUM+1)
    
        !///////�м������ˮλ��ֵ
        DO j=Element_NUM,1,-1
           !/////���ж������ٵ��ö�Ӧ�ĺ���ģ�飬�����Ԫ�����ˮ�ͬʱ�����Ԫ����С���������
           IF(Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
               Call CanalInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ELSEIF(Ele_Rel(j,2)==4)THEN  !/////������Ԫ��
               Call SideOutFlowInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ELSEIF(Ele_Rel(j,2)==6)THEN  !/////����բԪ��
               Call ControlGateInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ELSEIF(Ele_Rel(j,2)==7)THEN  !/////�����Ԫ��
               Call AreachangeInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ELSEIF(Ele_Rel(j,2)==10)THEN !/////������Ԫ��
               Call InvertedsiphonInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
	ELSEIF(Ele_Rel(j,2)==11)THEN !/////����Ԫ��
               Call BridgeInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ELSEIF(Ele_Rel(j,2)==12)THEN !/////��վԪ��
               Call PumpingstationInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ENDIF 
        END DO
 
    !***************************************!
    !***********�����ʼ������************!
    !***************************************!
    
    !//////����������///////////
    DO j=1,Element_NUM,1
        IF(Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
           L(j)=Ca_cl(Ele_Rel(j,3))%Incoord
           L(j+1)=Ca_cl(Ele_Rel(j,3))%Outcoord
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////������Ԫ��
           L(j)=Insi_cl(Ele_Rel(j,3))%Incoord
           L(j+1)=Insi_cl(Ele_Rel(j,3))%Outcoord
        ENDIF     
    END DO
    
    !//////�����������///////////
    DO j=1,Element_NUM,1
        IF(Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
           vc(j)=Q(j)/((Ca_cl(Ele_Rel(j,3))%BottomWidth+Ca_cl(Ele_Rel(j,3))%SideSlopes*h(j))*h(j))
           vc(j+1)=Q(j+1)/((Ca_cl(Ele_Rel(j,3))%BottomWidth+Ca_cl(Ele_Rel(j,3))%SideSlopes*h(j+1))*h(j+1))
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////������Ԫ��
           vc(j)=Q(j)/(Insi_cl(Ele_Rel(j,3))%ParallelNum*Insi_cl(Ele_Rel(j,3))%BottomWidth*h(j))
           vc(j+1)=Q(j+1)/(Insi_cl(Ele_Rel(j,3))%ParallelNum*Insi_cl(Ele_Rel(j,3))%BottomWidth*h(j+1))
        ENDIF     
    END DO
    
    !//////���������ˮλ��������õ���ˮ����ϸ�����ĵ׸߳�
    DO j=1,Element_NUM,1
       IF(Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
           Z(j)=Ca_cl(Ele_Rel(j,3))%InInvertElev+h(j)
           Z(j+1)=Ca_cl(Ele_Rel(j,3))%OutInvertElev+h(j+1)
       ELSEIF(Ele_Rel(j,2)==10)THEN !/////������Ԫ��
           Z(j)=Insi_cl(Ele_Rel(j,3))%InInvertElev+h(j)
           Z(j+1)=Insi_cl(Ele_Rel(j,3))%OutInvertElev+h(j+1)
       ENDIF 
    END DO
	Vol(Element_NUM+1)=0
    
    !/////�����ʼ����ˮ��������//////////
    OPEN(UNIT=20,FILE="OUTPUTFILE/InitialRESULT.txt")
        WRITE(20,"(A10,3X,A10,3X,A7,3X,A7,3X,A7,3X,A7,3X,A7,3X)")"�������","���","ˮ��","����","����","ˮλ","��ˮͷ"
        DO i=1,Element_NUM+1
            WRITE(20,"(I6,7X,F10.3,3X,F7.3,3X,F7.3,3X,F7.3,3X,F7.3,3X,F7.3,3X)")i,L(i),h(i),Q(i),vc(i),Z(i),Z(i)+vc(i)**2/2.0/GRAV
        END DO
    CLOSE(20)

    !******************************!
    !******�Ǻ㶨������Ԥ����******��
    !******************************!
    !����Ǻ㶨������ռ估����Ǻ㶨��������µ�Ԫ����Ϣ
    !//////����ϵ������ռ䣬֮�������ֵ���������������ģ��ĵ���
    !//////���㷽�̸�����ȡ����n(i)
    c_element=0
    DO i=1,Element_NUM,1
        c_element=c_element+nnn(i)
    END DO
    
    !//////////���µ�Ԫ�����б�ţ�����ȷ������Ԫ������ػ���������������������Ӱ�쵽�����
    !��Ҫ�Ƕ��������½��оֲ���ţ�����Ԫ�������������±��
    ALLOCATE(Ele_Rel_New(c_element,3))
    k=0
    gg=1
    DO i=1,Element_NUM
        DO j=1,nnn(i)
            k=k+1
            Ele_Rel_New(k,1)=k
            Ele_Rel_New(k,2)=Ele_Rel(i,2)
            IF(Ele_Rel(i,2)==5)THEN
                Ele_Rel_New(k,3)=j-1+gg
            ELSE
                Ele_Rel_New(k,3)=Ele_Rel(i,3)
            END IF
        END DO
        IF(Ele_Rel(i,2)==5)THEN
            gg=gg+nnn(i)
        END IF
    END DO
!    OPEN(UNIT=21,FILE="TEST.txt")
!        DO i=1,c_element
!            WRITE(21,*)(Ele_Rel_New(i,j),j=1,3)   
!        END DO
!    CLOSE(21)

    !ȷ���µļ���Ԫ���Ļ���������ʵ����ֻ��ȷ���µ���������Ϣ������Ԫ������ֱ�ӵȼ�
    Canal_New_NUM=gg-1
    ALLOCATE(n_canal(Canal_NUM))
    j=0
    DO i=1,Element_NUM
        IF(Ele_Rel(i,2)==5)THEN
            j=j+1
            n_canal(j)=nnn(i)
        END IF
    END DO
    ALLOCATE(Ca_cl_New(Canal_New_NUM))
    k=0
    DO i=1,Canal_NUM
        DO j=1,n_canal(i)
            k=k+1
            Ca_cl_New(k)%Serial_NUM=k
            Ca_cl_New(k)%Name=Ca_cl(i)%Name
            Ca_cl_New(k)%BottomWidth=Ca_cl(i)%BottomWidth
            Ca_cl_New(k)%SideSlopes=Ca_cl(i)%SideSlopes
            Ca_cl_New(k)%Incoord=Ca_cl(i)%Incoord+(j-1)*(Ca_cl(i)%Outcoord-Ca_cl(i)%Incoord)/n_canal(i)
            Ca_cl_New(k)%Outcoord=Ca_cl(i)%Incoord+j*(Ca_cl(i)%Outcoord-Ca_cl(i)%Incoord)/n_canal(i)
            Ca_cl_New(k)%InInvertElev=Ca_cl(i)%InInvertElev+(j-1)*(Ca_cl(i)%OutInvertElev-Ca_cl(i)%InInvertElev)/n_canal(i)
            Ca_cl_New(k)%OutInvertElev=Ca_cl(i)%InInvertElev+j*(Ca_cl(i)%OutInvertElev-Ca_cl(i)%InInvertElev)/n_canal(i)
            Ca_cl_New(k)%InTopElev=Ca_cl(i)%InTopElev+(j-1)*(Ca_cl(i)%OutTopElev-Ca_cl(i)%InTopElev)/n_canal(i)
            Ca_cl_New(k)%OutTopElev=Ca_cl(i)%InTopElev+j*(Ca_cl(i)%OutTopElev-Ca_cl(i)%InTopElev)/n_canal(i)
            Ca_cl_New(k)%Roughness=Ca_cl(i)%Roughness
        END DO
    END DO

    DEALLOCATE(h)
    DEALLOCATE(Q)
    DEALLOCATE(Z)
    DEALLOCATE(vc)
    DEALLOCATE(L)
End 
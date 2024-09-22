    SUBROUTINE Initial_Cal !wjb0413
    USE GLOBAL_VARIABLE

    IMPLICIT NONE

    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_main
    !TYPE(Config)::Config_main
    !TYPE(Initial_Result)::Initial_Result_main
    !TYPE(Condition)::Condition_0
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation�洢ϵͳ��š����ʹ��롢�ֲ����

    INTEGER::i,j,k
    INTEGER::ele_sg,DDM_sg    !��¼��Ԫ����λ�ã������Ŷգ�
    INTEGER::Element_NUM    !����Ԫ���ĸ���
    INTEGER,ALLOCATABLE::DDM_POSTION(:)  !����բ+��β����λ��
    REAL::Waterdepth_flag

    !/////�����Ϣ
    REAL,ALLOCATABLE::h(:),Q(:),Z(:),vc(:),L(:),Vol(:),Pool_Volume_ini(:)

    Element_NUM=Com_element_main.Element_NUM
    Ele_Rel=Com_element_main.Ele_Rel

    !/////��ʼ����
    !��������ˮλ����������洢�ռ�
    CALL Initial_Result_Auto_Update(Initial_Result_main,Com_element_main,"NUM")
    ALLOCATE(h(Element_NUM+1))
    ALLOCATE(Q(Element_NUM+1))
    ALLOCATE(Z(Element_NUM+1))
    ALLOCATE(vc(Element_NUM+1))
    ALLOCATE(L(Element_NUM+1))
    ALLOCATE(Vol(Element_NUM+1))
    ALLOCATE(Pool_Volume_ini(Com_element_main.Controlgate_NUM+1))
    !����Ԫ���洢��������洢�ռ�

    Waterdepth_flag=0

    !/////��ʼ��������������
    Q(1)=Boundary_cal.UpQ(1)
    j=0
    DO i=2,Element_NUM
        IF(Ele_Rel(i,2)==4)THEN
            j=j+1
            Q(i)=Q(i-1)-Boundary_cal.OutDischarge(j,1)
        Else
            Q(i)=Q(i-1)
        END IF
    END DO
    Q(Element_NUM+1)=Q(Element_NUM)

    !//////��ʼˮλ�ļ��㣬�Լ������εķֶΣ�Ŀ����ȷ���������ˮλ
    !//////ˮλ�ļ����Ǵ���������������
    !//////����ˮλ�߽����
    h(Element_NUM+1)=Boundary_cal.DownZ(1)-Com_element_main.Ca_cl(Com_element_main.Canal_NUM).OutInvertElev
    !h(Element_NUM+1)=Boundary_cal.DownZ(1)
    
    Waterdepth_flag=Waterdepth_flag+1
    Condition_0.h(Waterdepth_flag)=h(Element_NUM+1)
    Condition_0.Q(Waterdepth_flag)=Q(Element_NUM+1)

    !///////�м������ˮλ��ֵ
    DO j=Element_NUM,1,-1
        !/////���ж������ٵ��ö�Ӧ�ĺ���ģ�飬�����Ԫ�����ˮ�ͬʱ�����Ԫ����С���������
        IF(Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
            Call CanalInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==4)THEN  !/////������Ԫ��
            Call SideOutFlowInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==6)THEN  !/////����բԪ��
            Call ControlGateInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==7)THEN  !/////�����Ԫ��
            Call AreachangeInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////������Ԫ��
            Call InvertedsiphonInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==11)THEN !/////����Ԫ��
            Call BridgeInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==12)THEN !/////��վԪ��
            Call PumpingstationInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ENDIF
    END DO

    !***************************************!
    !***********�����ʼ������************!
    !***************************************!

    !//////����������///////////
    DO j=1,Element_NUM,1
        IF(Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
            L(j)=Com_element_main.Ca_cl(Ele_Rel(j,3)).Incoord
            L(j+1)=Com_element_main.Ca_cl(Ele_Rel(j,3)).Outcoord
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////������Ԫ��
            L(j)=Com_element_main.Insi_cl(Ele_Rel(j,3)).Incoord
            L(j+1)=Com_element_main.Insi_cl(Ele_Rel(j,3)).Outcoord
        ENDIF
    END DO

    !//////�����������///////////
    DO j=1,Element_NUM,1
        IF(Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
            vc(j)=Q(j)/((Com_element_main.Ca_cl(Ele_Rel(j,3)).BottomWidth+Com_element_main.Ca_cl(Ele_Rel(j,3)).SideSlopes*h(j))*h(j))
            vc(j+1)=Q(j+1)/((Com_element_main.Ca_cl(Ele_Rel(j,3)).BottomWidth+Com_element_main.Ca_cl(Ele_Rel(j,3)).SideSlopes*h(j+1))*h(j+1))
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////������Ԫ��
            vc(j)=Q(j)/(Com_element_main.Insi_cl(Ele_Rel(j,3)).ParallelNum*Com_element_main.Insi_cl(Ele_Rel(j,3)).BottomWidth*h(j))
            vc(j+1)=Q(j+1)/(Com_element_main.Insi_cl(Ele_Rel(j,3)).ParallelNum*Com_element_main.Insi_cl(Ele_Rel(j,3)).BottomWidth*h(j+1))
        ENDIF
    END DO

    !//////���������ˮλ��������õ���ˮ����ϸ�����ĵ׸߳�
    DO j=1,Element_NUM,1
        IF(Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
            Z(j)=Com_element_main.Ca_cl(Ele_Rel(j,3)).InInvertElev+h(j)
            Z(j+1)=Com_element_main.Ca_cl(Ele_Rel(j,3)).OutInvertElev+h(j+1)
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////������Ԫ��
            Z(j)=Com_element_main.Insi_cl(Ele_Rel(j,3)).InInvertElev+h(j)
            Z(j+1)=Com_element_main.Insi_cl(Ele_Rel(j,3)).OutInvertElev+h(j+1)
        ENDIF
    END DO
    Vol(Element_NUM+1)=0

    !///////ȷ��բ��+��λ�������ڵ�λ��//////
    ALLOCATE(DDM_POSTION(Com_element_main.Controlgate_NUM+2))
    DDM_sg=2
    DDM_POSTION(1)=1
    DO j=2,Element_NUM
        IF(Ele_Rel(j,2)==6)THEN
            DDM_POSTION(DDM_sg)=j
            DDM_sg=DDM_sg+1
        END IF
    END DO
    DDM_POSTION(Com_element_main.Controlgate_NUM+2)=Element_NUM

    !�����صĳ�ʼ�����ܺͼ���
    DO i=1,Com_element_main.Controlgate_NUM+1
        Pool_Volume_ini(i)=0.0
        DO j=DDM_POSTION(i),DDM_POSTION(i+1)-1
            Pool_Volume_ini(i)=Pool_Volume_ini(i)+Vol(j)
        END DO
    END DO

    Initial_Result_main.h=h
    Initial_Result_main.Q=Q
    Initial_Result_main.Z=Z
    Initial_Result_main.vc=vc
    Initial_Result_main.L=L
    Initial_Result_main.Vol=Vol
    Initial_Result_main.Pool_Volume=Pool_Volume_ini

    End
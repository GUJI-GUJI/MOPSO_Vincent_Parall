    SUBROUTINE Initial_Cal !wjb0413
    USE GLOBAL_VARIABLE

    IMPLICIT NONE

    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_main
    !TYPE(Config)::Config_main
    !TYPE(Initial_Result)::Initial_Result_main
    !TYPE(Condition)::Condition_0
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation存储系统编号、类型代码、局部编号

    INTEGER::i,j,k
    INTEGER::ele_sg,DDM_sg    !记录各元件的位置（除了桥墩）
    INTEGER::Element_NUM    !计算元件的个数
    INTEGER,ALLOCATABLE::DDM_POSTION(:)  !节制闸+首尾断面位置
    REAL::Waterdepth_flag

    !/////结果信息
    REAL,ALLOCATABLE::h(:),Q(:),Z(:),vc(:),L(:),Vol(:),Pool_Volume_ini(:)

    Element_NUM=Com_element_main.Element_NUM
    Ele_Rel=Com_element_main.Ele_Rel

    !/////初始计算
    !给各断面水位和流量分配存储空间
    CALL Initial_Result_Auto_Update(Initial_Result_main,Com_element_main,"NUM")
    ALLOCATE(h(Element_NUM+1))
    ALLOCATE(Q(Element_NUM+1))
    ALLOCATE(Z(Element_NUM+1))
    ALLOCATE(vc(Element_NUM+1))
    ALLOCATE(L(Element_NUM+1))
    ALLOCATE(Vol(Element_NUM+1))
    ALLOCATE(Pool_Volume_ini(Com_element_main.Controlgate_NUM+1))
    !给各元件存储个数分配存储空间

    Waterdepth_flag=0

    !/////初始各断面流量计算
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

    !//////初始水位的计算，以及各渠段的分段！目的是确定各断面的水位
    !//////水位的计算是从下游往上游推求
    !//////下游水位边界计算
    h(Element_NUM+1)=Boundary_cal.DownZ(1)-Com_element_main.Ca_cl(Com_element_main.Canal_NUM).OutInvertElev
    !h(Element_NUM+1)=Boundary_cal.DownZ(1)
    
    Waterdepth_flag=Waterdepth_flag+1
    Condition_0.h(Waterdepth_flag)=h(Element_NUM+1)
    Condition_0.Q(Waterdepth_flag)=Q(Element_NUM+1)

    !///////中间各断面水位赋值
    DO j=Element_NUM,1,-1
        !/////先判断类型再调用对应的函数模块，计算各元件入口水深，同时计算各元件的小断面个数。
        IF(Ele_Rel(j,2)==5)THEN  !/////渠道元件
            Call CanalInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==4)THEN  !/////分数口元件
            Call SideOutFlowInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==6)THEN  !/////节制闸元件
            Call ControlGateInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==7)THEN  !/////渐变段元件
            Call AreachangeInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
            Call InvertedsiphonInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==11)THEN !/////桥梁元件
            Call BridgeInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ELSEIF(Ele_Rel(j,2)==12)THEN !/////泵站元件
            Call PumpingstationInitialStageCalc(j,Waterdepth_flag,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),Vol(j))
        ENDIF
    END DO

    !***************************************!
    !***********输出初始计算结果************!
    !***************************************!

    !//////各断面的里程///////////
    DO j=1,Element_NUM,1
        IF(Ele_Rel(j,2)==5)THEN  !/////渠道元件
            L(j)=Com_element_main.Ca_cl(Ele_Rel(j,3)).Incoord
            L(j+1)=Com_element_main.Ca_cl(Ele_Rel(j,3)).Outcoord
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
            L(j)=Com_element_main.Insi_cl(Ele_Rel(j,3)).Incoord
            L(j+1)=Com_element_main.Insi_cl(Ele_Rel(j,3)).Outcoord
        ENDIF
    END DO

    !//////各断面的流速///////////
    DO j=1,Element_NUM,1
        IF(Ele_Rel(j,2)==5)THEN  !/////渠道元件
            vc(j)=Q(j)/((Com_element_main.Ca_cl(Ele_Rel(j,3)).BottomWidth+Com_element_main.Ca_cl(Ele_Rel(j,3)).SideSlopes*h(j))*h(j))
            vc(j+1)=Q(j+1)/((Com_element_main.Ca_cl(Ele_Rel(j,3)).BottomWidth+Com_element_main.Ca_cl(Ele_Rel(j,3)).SideSlopes*h(j+1))*h(j+1))
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
            vc(j)=Q(j)/(Com_element_main.Insi_cl(Ele_Rel(j,3)).ParallelNum*Com_element_main.Insi_cl(Ele_Rel(j,3)).BottomWidth*h(j))
            vc(j+1)=Q(j+1)/(Com_element_main.Insi_cl(Ele_Rel(j,3)).ParallelNum*Com_element_main.Insi_cl(Ele_Rel(j,3)).BottomWidth*h(j+1))
        ENDIF
    END DO

    !//////计算各断面水位，即计算得到的水深加上各断面的底高程
    DO j=1,Element_NUM,1
        IF(Ele_Rel(j,2)==5)THEN  !/////渠道元件
            Z(j)=Com_element_main.Ca_cl(Ele_Rel(j,3)).InInvertElev+h(j)
            Z(j+1)=Com_element_main.Ca_cl(Ele_Rel(j,3)).OutInvertElev+h(j+1)
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
            Z(j)=Com_element_main.Insi_cl(Ele_Rel(j,3)).InInvertElev+h(j)
            Z(j+1)=Com_element_main.Insi_cl(Ele_Rel(j,3)).OutInvertElev+h(j+1)
        ENDIF
    END DO
    Vol(Element_NUM+1)=0

    !///////确定闸门+首位断面所在的位置//////
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

    !各渠池的初始蓄量总和计算
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
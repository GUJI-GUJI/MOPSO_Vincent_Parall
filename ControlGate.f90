    !/////节制闸 元件代码 6
    !////////分成两个部分，一个是恒定状态的，一个是非恒定状态的
    !////初始系数
    SUBROUTINE ControlGateInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
    USE GLOBAL_VARIABLE

    IMPLICIT NONE
    !TYPE(Com_element)::Com_element_main
    !TYPE(Condition)::Condition_0
    REAL::Waterdepth_flag
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation存储系统编号、类型代码、局部编号
    !/////子程序与主调程序之间的数据交换
    INTEGER::j_NUM  !元件系统编号
    INTEGER::Serial !对应节制闸的局部编号
    REAL::h2,Q2,Q1,h1   !分别表示循环过程中进口断面和出口断面的水深和流量
    INTEGER::i
    INTEGER::nn
    Real::Forebay_h !前池水位和开度
    REAL::vvol  !元件蓄量

    Forebay_h=Com_element_main.Cogate_cl(Serial).Initial_gate_Zup-Com_element_main.Cogate_cl(Serial).InvertElev
    !///////节制闸分段数
    Com_element_main.n(j_NUM)=1!///永远只分一段
    h1=Forebay_h
    Waterdepth_flag=Waterdepth_flag+1
    Condition_0.h(Waterdepth_flag)=h1
    Condition_0.Q(Waterdepth_flag)=Q2
    vvol=h2*Com_element_main.Cogate_cl(Serial).Lenth*Com_element_main.Cogate_cl(Serial).RUNNum*Com_element_main.Cogate_cl(Serial).SingleWidth
    END


    !////////////////非恒定流系数
    Subroutine ControlGateDynamicCoeff(j_NUM,Serial,Opde,Qo2,Qo1,Zo2,Zo1,S2,T2,P2,V2,ele_vvol,Condition_cal)
    USE GLOBAL_VARIABLE

    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_cal
    TYPE(Condition)::Condition_cal
    INTEGER::j_NUM  !元件系统编号
    INTEGER::Serial !对应渠道的局部编号
    REAL::Qo2,Qo1,Zo2,Zo1 !分别表示n时刻j+1断面流量，j断面流量，j+1断面水位，j断面水位
    REAL::S2,T2,P2,V2   !追赶系数
    REAL::C,D,E,F,G,M   !线性方程组的系数
    REAL::Y1,Y2,Y3,Y4   !中间系数
    REAL::hh1,hh2   !水深
    REAL::BB  !过闸水流宽度
    REAL::Opde  !闸门开度
    REAL::CC    !综合系数
    REAL::ele_vvol


    !/////需要用到的水力要素预处理
    hh1=Zo1-Com_element_cal.Cogate_cl(Serial).InvertElev
    hh2=Zo2-Com_element_cal.Cogate_cl(Serial).InvertElev
    BB=Com_element_cal.Cogate_cl(Serial).RUNNum*Com_element_cal.Cogate_cl(Serial).SingleWidth
    ele_vvol=hh2*Com_element_cal.Cogate_cl(Serial).Lenth*Com_element_cal.Cogate_cl(Serial).RUNNum*Com_element_cal.Cogate_cl(Serial).SingleWidth
    if ((hh1<hh2)) then
        write(*,*)
        write(*,*)
        write(*,*)"=======WARNING（ControlGate.f90——非恒定流计算）====="
        write(*,*) "第",j_NUM,"个计算元件的"
        write(*,*) "节制闸闸前水位低于闸后水位"
        write(*,*) "说明上一时刻计算的水深不合理，尝试对其自动修正"
        hh2=hh1-0.05
    end if

    !!固定系数
    CC=(Com_element_cal.Cogate_cl(Serial).Discoeff)*BB*Opde*SQRT(2.0*9.81)   !模拟用 固定系数
    !!闸门开度-流量系数
    !CC=(Com_element_cal.Cogate_cl(Serial).G_A*Opde**2+Com_element_cal.Cogate_cl(Serial).G_B*Opde+Com_element_cal.Cogate_cl(Serial).G_C)*BB*Opde*SQRT(2.0*9.81)  !变系数
    !!闸门开度/水头差-流量系数
    !CC=(Com_element_cal.Cogate_cl(Serial).G_A*(Opde/(Zo1-Zo2))**2+Com_element_cal.Cogate_cl(Serial).G_B*(Opde/(Zo1-Zo2))+Com_element_cal.Cogate_cl(Serial).G_C)*BB*Opde*SQRT(2.0*9.81)




    !IF(Boundary_cal.Upsign==0)THEN
    !    !//////////中间系数
    !    Y1=Condition_cal.VV(j_NUM)+C
    !    Y2=F+E*Condition_cal.VV(j_NUM)
    !    Y3=D+Condition_cal.PP(j_NUM)
    !    Y4=M-E*Condition_cal.PP(j_NUM)
    !
    !    !////////追赶系数
    !    S2=(G*Y3-Y4)/(Y1*G+Y2)
    !    T2=(C*G-F)/(Y1*G+Y2)
    !    P2=Y3-Y1*S2
    !    V2=C-Y1*T2
    !ELSE
    !    !////////追赶系数
    !    S2=0.0
    !    T2=-1.0
    !    P2=M/F+Condition_cal.PP(j_NUM)
    !    V2=Condition_cal.VV(j_NUM)+1.0/F
    !ENDIF


    !追赶系数
    IF (Opde<=MAX(hh1,hh2))Then !闸门节制
        !避免闸前水位相差过小
        IF(abs(hh1-hh2)<0.01)Then
            hh1=hh1+0.02
            hh2=hh2-0.02
        ENDIF

        !避免闸前水位低于闸后水位
        IF(hh1<hh2)THEN
            hh1=hh2+0.05
            Write(*,'(*(G0.5))') "节制闸出现了闸前水位低于闸后水位的情况，已强行修正，请注意验证结果可靠性！"
        ENDIF

        !//////线性方程组系数赋值(仅考虑了淹没出流）
        C=0.0
        D=0.0
        E=0.0
        G=1.0
        F=CC/2.0/SQRT(hh1-hh2)
        M=CC*SQRT(hh1-hh2)/2.0
        IF(Boundary_Cal.Upsign==0)THEN
            Y1=Condition_cal.VV(j_NUM)+C
            Y2=F+E*Condition_cal.VV(j_NUM)
            Y3=D+Condition_cal.PP(j_NUM)
            Y4=M-E*Condition_cal.PP(j_NUM)

            S2=(G*Y3-Y4)/(Y1*G+Y2)
            T2=(C*G-F)/(Y1*G+Y2)
            P2=Y3-Y1*S2
            V2=C-Y1*T2
        ELSE
            S2=0.0
            T2=-1.0
            P2=M/F+Condition_cal.PP(j_NUM)
            V2=Condition_cal.VV(j_NUM)+1.0/F
        ENDIF


    ELSEIF(Opde>MAX(hh1,hh2))Then!闸门全开
        !Condition_cal.SS(j_NUM+1)=0.0
        !Condition_cal.TT(j_NUM+1)=-1.0
        !Condition_cal.PP(j_NUM+1)=Condition_cal.PP(j_NUM)
        !Condition_cal.VV(j_NUM+1)=Condition_cal.VV(j_NUM)
        S2=0
        T2=-1.0
        P2=Condition_cal.PP(j_NUM)
        V2=Condition_cal.VV(j_NUM)
    ENDIF






    END
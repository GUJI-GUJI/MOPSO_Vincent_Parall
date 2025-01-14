!/////桥梁 元件代码11
!////////分成两个部分，一个是恒定状态的，一个是非恒定状态的
!////初始系数
SUBROUTINE BridgeInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    
    !TYPE(Com_element)::Com_element_main
    !TYPE(Condition)::Condition_0
    REAL::Waterdepth_flag
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation存储系统编号、类型代码、局部编号
    
    !/////子程序与主调程序之间的数据交换
    INTEGER::j_NUM  !元件系统编号    
    INTEGER::Serial !对应倒虹吸的局部编号
    REAL::h2,Q2,Q1,h1   !分别表示循环过程中进口断面和出口断面的水深和流量
    INTEGER::nn
    !/////子程序内部计算涉及到的参数
    REAL::BB1,BB2   !过水断面宽度
    REAL::Ar_Brid   !桥墩处的面积
    REAL::Ar_NoBrid !桥墩后的面积
    REAL::u_Brid    !桥墩后的流速
    REAL::alpha     !A'/A 桥墩阻水比
    REAL::k         !桥墩形状系数
    REAL::vvol
    
    Ele_Rel=Com_element_main.Ele_Rel
    
    !///////桥梁分段数

    Com_element_main.n(j_NUM)=1
    k=0.9   !///桥墩形状系数

    !///////下游水深计算
    IF(Ele_Rel(j_NUM-1,2)==5.AND.Ele_Rel(j_NUM+1,2)==5)THEN
        h1=h2
        BB2=Com_element_main.Brid_cl(Serial).BottomWidth+2*Com_element_main.Brid_cl(Serial).SideSlopes*h2
        Ar_NoBrid=(BB2+Com_element_main.Ca_cl(Ele_Rel(j_NUM+1,3)).BottomWidth)*h2/2.0
        Ar_Brid=Com_element_main.Brid_cl(Serial).HinderWidth*h2
        u_Brid=Q2/Ar_NoBrid
        alpha=Ar_Brid/Ar_NoBrid
        h1=h2+2.0*k*(k+10*((u_Brid**2.0)/(2.0*9.81*h2))-0.6)*(alpha+15.0*alpha**4.0)*((u_Brid**2.0)/(2.0*9.81))
    END IF 
    
    Waterdepth_flag=Waterdepth_flag+1
    Condition_0.h(Waterdepth_flag)=h1
    Condition_0.Q(Waterdepth_flag)=Q2
    vvol=0
END


!////////////////非恒定流系数
Subroutine BridgeDynamicCoeff(j_NUM,Serial,Qo2,Qo1,Zo2,Zo1,S2,T2,P2,V2,ele_vvol,Condition_cal)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_cal
    TYPE(Condition)::Condition_cal

    INTEGER::j_NUM  !元件系统编号    
    INTEGER::Serial !对应倒虹吸的局部编号
    REAL::Qo2,Qo1,Zo2,Zo1 !分别表示n时刻j+1断面流量，j断面流量，j+1断面水位，j断面水位
    REAL::S2,T2,P2,V2   !追赶系数
    REAL::C,D,E,F,G,M   !线性方程组的系数
    REAL::Y1,Y2,Y3,Y4   !中间系数
    REAL::hh1,hh2   !水深
    REAL::BB1,BB2   !水面宽度
    REAL::Ar1,Ar2,Ar_Brid      !过水面积
    REAL::L_wet     !湿周
    REAL::R !水力半径
    REAL::k !桥墩形状系数
    REAL::alpha     !A'/A 桥墩阻水比
    REAL::u1,u2 !断面流速
    REAL::ele_vvol  !元件蓄量

    !/////需要用到的水力要素预处理
    k=0.9
    IF(Com_element_cal.Ele_Rel(j_NUM-1,2)==5.AND.Com_element_cal.Ele_Rel(j_NUM+1,2)==5)THEN
        hh1=Zo1-Com_element_cal.Brid_cl(Serial).InInvertElev
        hh2=Zo2-Com_element_cal.Brid_cl(Serial).OutInvertElev
        BB1=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).BottomWidth+2*Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).SideSlopes*hh1
        BB2=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).BottomWidth+2*Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).SideSlopes*hh2
        Ar1=(BB1+Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).BottomWidth)*hh1/2.0
        Ar2=(BB2+Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).BottomWidth)*hh2/2.0
        Ar_Brid=Com_element_cal.Brid_cl(Serial).HinderWidth*hh2
        alpha=Ar_Brid/Ar2
        ele_vvol=0
    
        u1=Qo1/Ar1
        u2=Qo2/Ar2
    END IF

    !//////线性方程组系数赋值
    C=0.0
    D=0.0
    E=0.0
    G=0.0
    F=1.0
    M=-2*k*(k+10*((u2**2.0)/(2.0*9.81*hh2))-0.6)*(alpha+15.0*alpha**4.0)*((u2**2.0)/(2.0*9.81))

    
    IF(Boundary_cal.Upsign==0)THEN
    !//////////中间系数
        Y1=Condition_cal.VV(j_NUM)+C
        Y2=F+E*Condition_cal.VV(j_NUM)
        Y3=D+Condition_cal.PP(j_NUM)
        Y4=M-E*Condition_cal.PP(j_NUM)

        !////////追赶系数
        S2=(G*Y3-Y4)/(Y1*G+Y2)
        T2=(C*G-F)/(Y1*G+Y2)
        P2=Y3-Y1*S2
        V2=C-Y1*T2
    ELSE
        !////////追赶系数
        S2=0.0
        T2=-1.0
        P2=Condition_cal.PP(j_NUM)+M
        V2=Condition_cal.VV(j_NUM)
    ENDIF
END 
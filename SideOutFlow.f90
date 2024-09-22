!/////分水口 元件代码 4
!////////分成两个部分，一个是恒定状态的，一个是非恒定状态的
!////初始系数
SUBROUTINE SideOutFlowInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    
    !TYPE(Com_element)::Com_element_main
    !TYPE(Condition)::Condition_0
    REAL::Waterdepth_flag
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation存储系统编号、类型代码、局部编号
    !/////子程序与主调程序之间的数据交换
    INTEGER::j_NUM  !元件系统编号    
    INTEGER::Serial !对应分水口的局部编号
    REAL::h2,Q2,Q1,h1   !分别表示循环过程中进口断面和出口断面的水深和流量
    INTEGER::nn
    INTEGER::i,k
    REAL::vvol  !元件蓄量
    
    !/////子程序内部计算涉及到的参数
    REAL::BB1,BB2   !过水断面宽度
    REAL::Ar1,Ar2   !过水面积
    REAL::u1,u2 !断面流速
    REAL::aa,ee,dhh_S
    
    Ele_Rel=Com_element_main.Ele_Rel
    
    !///////分水口分段数
    
    Com_element_main.n(j_NUM)=1!///永远只分一段
    !///////下游水深计算
    h1=h2
    Do k=1,200,1
        aa=h1
        BB1=Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).BottomWidth+2*Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).SideSlopes*h1
        BB2=Com_element_main.Ca_cl(Ele_Rel(j_NUM+1,3)).BottomWidth+2*Com_element_main.Ca_cl(Ele_Rel(j_NUM+1,3)).SideSlopes*h2
        Ar1=(BB1+Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).BottomWidth)*h1/2.0
        Ar2=(BB2+Com_element_main.Ca_cl(Ele_Rel(j_NUM+1,3)).BottomWidth)*h2/2.0
        u1=Q1/Ar1
        u2=Q2/Ar2
        h1=h2+(u2**2-u1**2)/2.0/9.81 !流速水头差
        IF(ABS(h1-aa)<2.22*2.718281828459**(-16))THEN   !在满足假设和实际相差不大的情况下，进行判断
            EXIT
        END IF  
    END DO

    Waterdepth_flag=Waterdepth_flag+1
    Condition_0.h(Waterdepth_flag)=h1
    Condition_0.Q(Waterdepth_flag)=Q2
    vvol=(h1+h2)/2.0*((Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).BottomWidth+Com_element_main.Ca_cl(Ele_Rel(j_NUM+1,3)).BottomWidth)/2.0)*10.0 &         !固定分水口长度为10m
         +Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).SideSlopes*h1*h2*10.0 &
         +1.0/3.0*(Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).SideSlopes*(h1-h2)**2)*10.0
END

!////////////////非恒定流系数
!////////////////非恒定流系数
Subroutine SideOutFlowDynamicCoeff(j_NUM,Serial,Qout,S2,T2,P2,V2,ele_vvol,Condition_cal)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_cal
    TYPE(Condition)::Condition_cal
    INTEGER::j_NUM  !元件系统编号    
    INTEGER::Serial !对应渠道的局部编号
    REAL::Qout  !分水流量
    REAL::S2,T2,P2,V2   !追赶系数
    REAL::C,D,E,F,G,M   !线性方程组的系数
    REAL::Y1,Y2,Y3,Y4   !中间系数
    REAL::ele_vvol
    
    ele_vvol=0

    !/////需要用到的水力要素预处理
    !//////线性方程组系数赋值
    C=0.0
    D=-Qout
    E=0.0
    G=0.0
    F=1.0
    M=0.0
    
    IF(Boundary_cal.Upsign==0)THEN
        !//////////中间系数
        Y1=Condition_cal.VV(j_NUM)+C
        Y2=F+E*Condition_cal.PP(j_NUM)
        Y3=D+Condition_cal.PP(j_NUM)
        Y4=M-E*Condition_cal.PP(j_NUM)

        !////////追赶系数
        S2=(G*Y3-Y4)/(Y1*G+Y2)
        T2=(C*G-F)/(Y1*G+Y2)
        P2=Y3-Y1*S2
        V2=C-Y1*T2
    ELSE
        !////////追赶系数
        S2=Qout
        T2=-1.0
        P2=Condition_cal.PP(j_NUM)+Condition_cal.VV(j_NUM)*(-Qout)
        V2=Condition_cal.VV(j_NUM)
    END IF
END 
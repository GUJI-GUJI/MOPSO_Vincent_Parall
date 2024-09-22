!/////倒虹吸 元件代码 10
!////////分成两个部分，一个是恒定状态的，一个是非恒定状态的
!////初始系数
SUBROUTINE InvertedsiphonInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
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
    INTEGER::i
    INTEGER::nn
    REAL::vvol
    !/////子程序内部计算涉及到的参数
    REAL::BB1,BB2   !过水断面宽度
    REAL::Ar_Insi   !虹吸管的面积
    REAL::u_Insi    !虹吸管的流速
    REAL::L_wet !虹吸管湿周
    REAL::Hydr_R    !水力半径 !
    REAL::K_Insi    !流量模数
    
    Ele_Rel=Com_element_main.Ele_Rel
    
    !///////倒虹吸分段数  
    Com_element_main.n(j_NUM)=1!///永远只分一段

    !///////下游水深计算
    h1=h2
    Ar_Insi=Com_element_main.Insi_cl(Serial).ParallelNum*Com_element_main.Insi_cl(Serial).BottomWidth*Com_element_main.Insi_cl(Serial).Height
    BB1=Com_element_main.Insi_cl(Serial).ParallelNum*Com_element_main.Insi_cl(Serial).BottomWidth
    BB2=Com_element_main.Insi_cl(Serial).ParallelNum*Com_element_main.Insi_cl(Serial).BottomWidth
    L_wet=Com_element_main.Insi_cl(Serial).ParallelNum*(Com_element_main.Insi_cl(Serial).BottomWidth*2+2*Com_element_main.Insi_cl(Serial).Height)
    Hydr_R=Ar_Insi/L_wet
    u_Insi=Q2/Ar_Insi
    K_Insi=1.0/(Com_element_main.Insi_cl(Serial).Roughness)*Hydr_R**(1.0/6.0)*Ar_Insi*Hydr_R**(1.0/6.0)
    
    !利用能量守恒方程计算水深后计算水位
    h1=Com_element_main.Insi_cl(Serial).OutInvertElev+h2-Com_element_main.Insi_cl(Serial).InInvertElev&
        &+Com_element_main.Insi_cl(Serial).Insi_coeff*Com_element_main.Insi_cl(Serial).Length/4.0/Hydr_R*u_Insi**2.0/2.0/9.81

    Waterdepth_flag=Waterdepth_flag+1
    Condition_0.h(Waterdepth_flag)=h1
    Condition_0.Q(Waterdepth_flag)=Q2
    vvol=(Com_element_main.Insi_cl(Serial).BottomWidth*Com_element_main.Insi_cl(Serial).ParallelNum-(Com_element_main.Insi_cl(Serial).ParallelNum-1)*&
        Com_element_main.Insi_cl(Serial).BottomWidth)*Com_element_main.Insi_cl(Serial).Height*Com_element_main.Insi_cl(Serial).Length
END

!////////////////非恒定流系数
!////////////////非恒定流系数
Subroutine InvertedsiphonDynamicCoeff(j_NUM,Serial,Qo2,Qo1,Zo2,Zo1,S2,T2,P2,V2,ele_vvol,Condition_cal)
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
    REAL::Ar1,Ar2,Ar_Insi      !过水面积
    REAL::L_wet     !湿周
    REAL::R !水力半径
    REAL::K !流量模数
    REAL::C1,C2 !谢才系数
    REAL::u1,u2,u_Insi !断面流速
    REAL::ele_vvol

    !/////需要用到的水力要素预处理
    hh1=Zo1-Com_element_cal.Insi_cl(Serial).InInvertElev
    hh2=Zo2-Com_element_cal.Insi_cl(Serial).OutInvertElev
    BB1=Com_element_cal.Insi_cl(Serial).ParallelNum*Com_element_cal.Insi_cl(Serial).BottomWidth
    BB2=Com_element_cal.Insi_cl(Serial).ParallelNum*Com_element_cal.Insi_cl(Serial).BottomWidth
    Ar1=BB1*hh1
    Ar2=BB2*hh2
    Ar_Insi=Com_element_cal.Insi_cl(Serial).ParallelNum*Com_element_cal.Insi_cl(Serial).BottomWidth*Com_element_cal.Insi_cl(Serial).Height
    L_wet=Com_element_cal.Insi_cl(Serial).ParallelNum*2.0*(Com_element_cal.Insi_cl(Serial).BottomWidth+Com_element_cal.Insi_cl(Serial).Height)
    ele_vvol=(Com_element_cal.Insi_cl(Serial).BottomWidth*Com_element_cal.Insi_cl(Serial).ParallelNum-(Com_element_cal.Insi_cl(Serial).ParallelNum-1)*&
        Com_element_cal.Insi_cl(Serial).BottomWidth)*Com_element_cal.Insi_cl(Serial).Height*Com_element_cal.Insi_cl(Serial).Length

    u1=Qo1/Ar1
    u2=Qo2/Ar2
    u_Insi=(Qo1+Qo2)/2.0/Ar_Insi
    R=Ar_Insi/L_wet
    K=1.0/(Com_element_cal.Insi_cl(Serial).Roughness)*R**(1.0/6.0)*Ar_Insi*R**(1.0/6.0)
    !//////线性方程组系数赋值
    C=0.0
    D=0.0
    E=0.0
    G=0.0
    F=1.0

    !M=-(Insi_cl(Serial).Insi_coeff*u_Insi**2.0/2.0/GRAV+Qo1**2.0/K**2.0*Insi_cl(Serial).Length)
    M=-(Com_element_cal.Insi_cl(Serial).Insi_coeff*Com_element_cal.Insi_cl(Serial).Length/4.0/R*u_Insi**2.0/2.0/9.81)
    !M=-(Insi_cl(Serial).Insi_coeff)
    !M的表达式，也就是倒虹吸的水里损失到底是多少？怎么表示，系数是多少，是一个值得研究的问题。
    
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
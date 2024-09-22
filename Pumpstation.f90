!/////泵站 元件代码 12
!////////分成两个部分，一个是恒定状态的，一个是非恒定状态的

!////初始系数
SUBROUTINE PumpingstationInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Com_element)::Com_element_main
    !TYPE(Condition)::Condition_0
    REAL::Waterdepth_flag
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation存储系统编号、类型代码、局部编号
    !/////子程序与主调程序之间的数据交换
    INTEGER::j_NUM  !元件系统编号    
    INTEGER::Serial !对应泵站的局部编号
    REAL::h2,Q2,Q1,h1   !分别表示循环过程中进口断面和出口断面的水深和流量
    INTEGER::i
    INTEGER::nn
    Real::Forebay_h !前池水位和开度
    Real::vvol
    
    Ele_Rel=Com_element_main.Ele_Rel
    Forebay_h=PumpST_cl(Serial).Initial_pump_Zup-PumpST_cl(Serial).InvertElev
    !///////泵站分段数 
    Com_element_main.n(j_NUM)=1!///永远只分一段
    h1=Forebay_h  
    WRITE(23,"(F7.3,3X,F7.3,3X)")h1,Q2
    vvol=0
END

!////////////////非恒定流系数
Subroutine PumpingstationDynamicCoeff(j_NUM,Serial,Bladeangle_RT,length,Qo2,Qo1,Zo2,Zo1,S2,T2,P2,V2,ele_vvol,Condition_cal)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_cal
    TYPE(Condition)::Condition_cal
    INTEGER::j_NUM  !元件系统编号    
    INTEGER::Serial !对应泵站的局部编号
    REAL::Qo2,Qo1,Zo2,Zo1 !分别表示n时刻j+1断面流量，j断面流量，j+1断面水位，j断面水位
    REAL::S2,T2,P2,V2   !追赶系数
    REAL::C,D,E,F,G,M   !线性方程组的系数
    REAL::Y1,Y2,Y3,Y4   !中间系数
    REAL::hh1,hh2   !水深
    INTEGER::length
    REAL::Bladeangle_RT(length)
    REAL:: A_ps,B_ps,C_ps      !泵站组合特性曲线系数
    INTEGER::i  !转角序列号
    REAL:: ele_vvol
    
    
    !由水泵转角获取泵站曲线
    CALL PumpToPumpingstation(Bladeangle_RT,Serial,A_ps,B_ps,C_ps)
    
    !///////调用泵站组合转角系数的查询函数(先查询转角，再调用系数)
    !CALL PumpToPumpingstation(
    !DO i=1,PumpST_cl(Serial).Change_NUM-1
        !IF((Cal_Time-1)*DeltaT>=PumpST_cl(Serial).Bladeangle(i,1).AND.(Cal_Time-1)*DeltaT<PumpST_cl(Serial).Bladeangle(i+1,1))THEN
            
        !ELSEIF((Cal_Time-1)*DeltaT==PumpST_cl(Serial).Bladeangle(PumpST_cl(Serial).Change_NUM,1))THEN
            
        !ENDIF
    !END DO
    ele_vvol=0

    !//////线性方程组系数赋值
    C=0.0
    D=0.0
    E=0.0
    G=-(2*A_ps*Qo2+B_ps)
    F=1.0
    M=-A_ps*Qo2**2+C_ps
    
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
        P2=Condition_cal.PP(j_NUM)-A_ps*Qo2**2+C_ps
        V2=(2*A_ps*Qo2+B_ps)-Condition_cal.VV(j_NUM)
    ENDIF
END 
    
    
SUBROUTINE PumpToPumpingstation(Bladeangle_RT,Serial,A_ps,B_ps,C_ps)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !输入输出
    INTEGER,INTENT(IN)::Serial !泵站编号
    REAL,INTENT(IN)::Bladeangle_RT(PumpST_cl(Serial).UNIT_NUM) !输入的转角序列
    REAL,INTENT(OUT)::A_ps!输出的三个系数
    REAL,INTENT(OUT)::B_ps
    REAL,INTENT(OUT)::C_ps
    !计算中间值
    INTEGER,ALLOCATABLE::Bladeangle_Serial(:)!转角对应的序号
    INTEGER::PS_Bladeangle_Serial!转角组合在泵站中的序号
    REAL::Q_sum(100)!泵站的离散曲线
    REAL,ALLOCATABLE::Q_pump(:,:)!水泵的离散曲线
    REAL::H_Pump,H_Pump_min!取扬程范围最小的泵作为泵站扬程范围
    REAL::H_step!扬程离散步长
    REAL::Bladeangle_D,Bladeangle_D_min!当前转角与有曲线的转角的差，用于寻找当前转角对应曲线
    !标记及工具
    INTEGER::i,j,flag,weight
    REAL::a,b,c!为方便写公式，对应水泵的abc系数
    REAL::X_co(3,3)!系数矩阵，
    REAL::Y_co(3)!值矩阵
    REAL::A_co(3)!解矩阵
    DOUBLE PRECISION::X_co_n(3,3)!系数矩阵的逆矩阵
    REAL::k!求逆矩阵中间量
    
    !寻找当前转角对应曲线序号
    ALLOCATE(Bladeangle_Serial(PumpST_cl(Serial).UNIT_NUM))
    ALLOCATE(Q_pump(PumpST_cl(Serial).UNIT_NUM,100))
    DO i=1,PumpST_cl(Serial).UNIT_NUM
        Bladeangle_D_min=ABS(Bladeangle_RT(i)-Pump_cl(PumpST_cl(Serial).Pump_Type(i)).Bladeangle_Serial(1))
        Bladeangle_Serial(i)=1
        DO j=1,Pump_cl(PumpST_cl(Serial).Pump_Type(i)).Coefficient_NUM
            Bladeangle_D=ABS(Bladeangle_RT(i)-Pump_cl(PumpST_cl(Serial).Pump_Type(i)).Bladeangle_Serial(j))
            IF(Bladeangle_D<Bladeangle_D_min)THEN
                Bladeangle_D_min=Bladeangle_D
                Bladeangle_Serial(i)=j
            END IF
        END DO
    END DO
    !计算该转角组合在泵站中的序号
    PS_Bladeangle_Serial=0
    DO i=1,PumpST_cl(Serial).UNIT_NUM-1
        weight=1
        DO j=i+1,PumpST_cl(Serial).UNIT_NUM
            weight=weight*Pump_cl(PumpST_cl(Serial).Pump_Type(j)).Coefficient_NUM
        END DO
        PS_Bladeangle_Serial=PS_Bladeangle_Serial+Bladeangle_Serial(i)*weight
    END DO
    PS_Bladeangle_Serial=PS_Bladeangle_Serial+Bladeangle_Serial(PumpST_cl(Serial).UNIT_NUM)
    !如果已经计算完成，则跳出
    IF(ISNAN(PumpST_cl(Serial).Coefficient_A(PS_Bladeangle_Serial))==.false.)THEN
        A_ps=PumpST_cl(Serial).Coefficient_A(PS_Bladeangle_Serial)
        B_ps=PumpST_cl(Serial).Coefficient_B(PS_Bladeangle_Serial)
        C_ps=PumpST_cl(Serial).Coefficient_C(PS_Bladeangle_Serial)
        RETURN
    END IF
    !计算泵站扬程范围
    H_Pump_min=Pump_cl(PumpST_cl(Serial).Pump_Type(1)).Coefficient_C(Bladeangle_Serial(1))
    DO i=1,PumpST_cl(Serial).UNIT_NUM
        H_Pump=Pump_cl(PumpST_cl(Serial).Pump_Type(i)).Coefficient_C(Bladeangle_Serial(i))
        IF(H_Pump<H_Pump_min)THEN
            H_Pump_min=H_Pump
        END IF
    END DO
    !计算扬程离散步长
    H_step=H_Pump_min/100
    
    !离散各水泵曲线并相加
    Q_sum=0
    DO i=1,100
        DO j=1,PumpST_cl(Serial).UNIT_NUM
            a=Pump_cl(PumpST_cl(Serial).Pump_Type(j)).Coefficient_A(Bladeangle_Serial(j))
            b=Pump_cl(PumpST_cl(Serial).Pump_Type(j)).Coefficient_B(Bladeangle_Serial(j))
            c=Pump_cl(PumpST_cl(Serial).Pump_Type(j)).Coefficient_C(Bladeangle_Serial(j))-i*H_step!每次计算将曲线向下平移求根得到流量
            Q_pump(j,i)=(-b-sqrt(b*b-4*a*c))/(2*a)
            Q_sum(i)=Q_sum(i)+Q_pump(j,i)
        END DO
    END DO
    
    !!!!!!!!!!拟合曲线得到系数!!!!!!!!!!!!!
    !计算系数矩阵
    X_co=0
    DO i=1,100
        X_co(1,1)=X_co(1,1)+1
        X_co(1,2)=X_co(1,2)+Q_sum(i)
        X_co(1,3)=X_co(1,3)+Q_sum(i)**2
        X_co(2,2)=X_co(2,2)+Q_sum(i)**2
        X_co(2,3)=X_co(2,3)+Q_sum(i)**3
        X_co(3,3)=X_co(3,3)+Q_sum(i)**4
    END DO
    X_co(2,1)=X_co(1,2)
    X_co(3,1)=X_co(1,3)
    X_co(3,2)=X_co(2,3)
    !计算值矩阵
    Y_co=0
    DO i=1,100
        Y_co(1)=Y_co(1)+i*H_step
        Y_co(2)=Y_co(2)+Q_sum(i)*i*H_step
        Y_co(3)=Y_co(3)+i*H_step*(Q_sum(i)**2)
    END DO
    !系数矩阵求逆
    X_co_n=X_co
    CALL BRINV(X_co_n,3,flag)
    !求解
    A_co=matmul(X_co_n,Y_co)
    A_ps=A_co(3)
    B_ps=A_co(2)
    C_ps=A_co(1)
    PumpST_cl(Serial).Coefficient_A(PS_Bladeangle_Serial)=A_ps
    PumpST_cl(Serial).Coefficient_B(PS_Bladeangle_Serial)=B_ps
    PumpST_cl(Serial).Coefficient_C(PS_Bladeangle_Serial)=C_ps
END SUBROUTINE PumpToPumpingstation
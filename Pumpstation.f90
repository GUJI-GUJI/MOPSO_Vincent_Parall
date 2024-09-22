!/////��վ Ԫ������ 12
!////////�ֳ��������֣�һ���Ǻ㶨״̬�ģ�һ���ǷǺ㶨״̬��

!////��ʼϵ��
SUBROUTINE PumpingstationInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Com_element)::Com_element_main
    !TYPE(Condition)::Condition_0
    REAL::Waterdepth_flag
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation�洢ϵͳ��š����ʹ��롢�ֲ����
    !/////�ӳ�������������֮������ݽ���
    INTEGER::j_NUM  !Ԫ��ϵͳ���    
    INTEGER::Serial !��Ӧ��վ�ľֲ����
    REAL::h2,Q2,Q1,h1   !�ֱ��ʾѭ�������н��ڶ���ͳ��ڶ����ˮ�������
    INTEGER::i
    INTEGER::nn
    Real::Forebay_h !ǰ��ˮλ�Ϳ���
    Real::vvol
    
    Ele_Rel=Com_element_main.Ele_Rel
    Forebay_h=PumpST_cl(Serial).Initial_pump_Zup-PumpST_cl(Serial).InvertElev
    !///////��վ�ֶ��� 
    Com_element_main.n(j_NUM)=1!///��Զֻ��һ��
    h1=Forebay_h  
    WRITE(23,"(F7.3,3X,F7.3,3X)")h1,Q2
    vvol=0
END

!////////////////�Ǻ㶨��ϵ��
Subroutine PumpingstationDynamicCoeff(j_NUM,Serial,Bladeangle_RT,length,Qo2,Qo1,Zo2,Zo1,S2,T2,P2,V2,ele_vvol,Condition_cal)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_cal
    TYPE(Condition)::Condition_cal
    INTEGER::j_NUM  !Ԫ��ϵͳ���    
    INTEGER::Serial !��Ӧ��վ�ľֲ����
    REAL::Qo2,Qo1,Zo2,Zo1 !�ֱ��ʾnʱ��j+1����������j����������j+1����ˮλ��j����ˮλ
    REAL::S2,T2,P2,V2   !׷��ϵ��
    REAL::C,D,E,F,G,M   !���Է������ϵ��
    REAL::Y1,Y2,Y3,Y4   !�м�ϵ��
    REAL::hh1,hh2   !ˮ��
    INTEGER::length
    REAL::Bladeangle_RT(length)
    REAL:: A_ps,B_ps,C_ps      !��վ�����������ϵ��
    INTEGER::i  !ת�����к�
    REAL:: ele_vvol
    
    
    !��ˮ��ת�ǻ�ȡ��վ����
    CALL PumpToPumpingstation(Bladeangle_RT,Serial,A_ps,B_ps,C_ps)
    
    !///////���ñ�վ���ת��ϵ���Ĳ�ѯ����(�Ȳ�ѯת�ǣ��ٵ���ϵ��)
    !CALL PumpToPumpingstation(
    !DO i=1,PumpST_cl(Serial).Change_NUM-1
        !IF((Cal_Time-1)*DeltaT>=PumpST_cl(Serial).Bladeangle(i,1).AND.(Cal_Time-1)*DeltaT<PumpST_cl(Serial).Bladeangle(i+1,1))THEN
            
        !ELSEIF((Cal_Time-1)*DeltaT==PumpST_cl(Serial).Bladeangle(PumpST_cl(Serial).Change_NUM,1))THEN
            
        !ENDIF
    !END DO
    ele_vvol=0

    !//////���Է�����ϵ����ֵ
    C=0.0
    D=0.0
    E=0.0
    G=-(2*A_ps*Qo2+B_ps)
    F=1.0
    M=-A_ps*Qo2**2+C_ps
    
    IF(Boundary_cal.Upsign==0)THEN
        !//////////�м�ϵ��
        Y1=Condition_cal.VV(j_NUM)+C
        Y2=F+E*Condition_cal.VV(j_NUM)
        Y3=D+Condition_cal.PP(j_NUM)
        Y4=M-E*Condition_cal.PP(j_NUM)

        !////////׷��ϵ��
        S2=(G*Y3-Y4)/(Y1*G+Y2)
        T2=(C*G-F)/(Y1*G+Y2)
        P2=Y3-Y1*S2
        V2=C-Y1*T2
    ELSE
        !////////׷��ϵ��
        S2=0.0
        T2=-1.0
        P2=Condition_cal.PP(j_NUM)-A_ps*Qo2**2+C_ps
        V2=(2*A_ps*Qo2+B_ps)-Condition_cal.VV(j_NUM)
    ENDIF
END 
    
    
SUBROUTINE PumpToPumpingstation(Bladeangle_RT,Serial,A_ps,B_ps,C_ps)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !�������
    INTEGER,INTENT(IN)::Serial !��վ���
    REAL,INTENT(IN)::Bladeangle_RT(PumpST_cl(Serial).UNIT_NUM) !�����ת������
    REAL,INTENT(OUT)::A_ps!���������ϵ��
    REAL,INTENT(OUT)::B_ps
    REAL,INTENT(OUT)::C_ps
    !�����м�ֵ
    INTEGER,ALLOCATABLE::Bladeangle_Serial(:)!ת�Ƕ�Ӧ�����
    INTEGER::PS_Bladeangle_Serial!ת������ڱ�վ�е����
    REAL::Q_sum(100)!��վ����ɢ����
    REAL,ALLOCATABLE::Q_pump(:,:)!ˮ�õ���ɢ����
    REAL::H_Pump,H_Pump_min!ȡ��̷�Χ��С�ı���Ϊ��վ��̷�Χ
    REAL::H_step!�����ɢ����
    REAL::Bladeangle_D,Bladeangle_D_min!��ǰת���������ߵ�ת�ǵĲ����Ѱ�ҵ�ǰת�Ƕ�Ӧ����
    !��Ǽ�����
    INTEGER::i,j,flag,weight
    REAL::a,b,c!Ϊ����д��ʽ����Ӧˮ�õ�abcϵ��
    REAL::X_co(3,3)!ϵ������
    REAL::Y_co(3)!ֵ����
    REAL::A_co(3)!�����
    DOUBLE PRECISION::X_co_n(3,3)!ϵ������������
    REAL::k!��������м���
    
    !Ѱ�ҵ�ǰת�Ƕ�Ӧ�������
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
    !�����ת������ڱ�վ�е����
    PS_Bladeangle_Serial=0
    DO i=1,PumpST_cl(Serial).UNIT_NUM-1
        weight=1
        DO j=i+1,PumpST_cl(Serial).UNIT_NUM
            weight=weight*Pump_cl(PumpST_cl(Serial).Pump_Type(j)).Coefficient_NUM
        END DO
        PS_Bladeangle_Serial=PS_Bladeangle_Serial+Bladeangle_Serial(i)*weight
    END DO
    PS_Bladeangle_Serial=PS_Bladeangle_Serial+Bladeangle_Serial(PumpST_cl(Serial).UNIT_NUM)
    !����Ѿ�������ɣ�������
    IF(ISNAN(PumpST_cl(Serial).Coefficient_A(PS_Bladeangle_Serial))==.false.)THEN
        A_ps=PumpST_cl(Serial).Coefficient_A(PS_Bladeangle_Serial)
        B_ps=PumpST_cl(Serial).Coefficient_B(PS_Bladeangle_Serial)
        C_ps=PumpST_cl(Serial).Coefficient_C(PS_Bladeangle_Serial)
        RETURN
    END IF
    !�����վ��̷�Χ
    H_Pump_min=Pump_cl(PumpST_cl(Serial).Pump_Type(1)).Coefficient_C(Bladeangle_Serial(1))
    DO i=1,PumpST_cl(Serial).UNIT_NUM
        H_Pump=Pump_cl(PumpST_cl(Serial).Pump_Type(i)).Coefficient_C(Bladeangle_Serial(i))
        IF(H_Pump<H_Pump_min)THEN
            H_Pump_min=H_Pump
        END IF
    END DO
    !���������ɢ����
    H_step=H_Pump_min/100
    
    !��ɢ��ˮ�����߲����
    Q_sum=0
    DO i=1,100
        DO j=1,PumpST_cl(Serial).UNIT_NUM
            a=Pump_cl(PumpST_cl(Serial).Pump_Type(j)).Coefficient_A(Bladeangle_Serial(j))
            b=Pump_cl(PumpST_cl(Serial).Pump_Type(j)).Coefficient_B(Bladeangle_Serial(j))
            c=Pump_cl(PumpST_cl(Serial).Pump_Type(j)).Coefficient_C(Bladeangle_Serial(j))-i*H_step!ÿ�μ��㽫��������ƽ������õ�����
            Q_pump(j,i)=(-b-sqrt(b*b-4*a*c))/(2*a)
            Q_sum(i)=Q_sum(i)+Q_pump(j,i)
        END DO
    END DO
    
    !!!!!!!!!!������ߵõ�ϵ��!!!!!!!!!!!!!
    !����ϵ������
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
    !����ֵ����
    Y_co=0
    DO i=1,100
        Y_co(1)=Y_co(1)+i*H_step
        Y_co(2)=Y_co(2)+Q_sum(i)*i*H_step
        Y_co(3)=Y_co(3)+i*H_step*(Q_sum(i)**2)
    END DO
    !ϵ����������
    X_co_n=X_co
    CALL BRINV(X_co_n,3,flag)
    !���
    A_co=matmul(X_co_n,Y_co)
    A_ps=A_co(3)
    B_ps=A_co(2)
    C_ps=A_co(1)
    PumpST_cl(Serial).Coefficient_A(PS_Bladeangle_Serial)=A_ps
    PumpST_cl(Serial).Coefficient_B(PS_Bladeangle_Serial)=B_ps
    PumpST_cl(Serial).Coefficient_C(PS_Bladeangle_Serial)=C_ps
END SUBROUTINE PumpToPumpingstation
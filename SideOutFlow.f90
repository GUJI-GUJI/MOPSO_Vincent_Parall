!/////��ˮ�� Ԫ������ 4
!////////�ֳ��������֣�һ���Ǻ㶨״̬�ģ�һ���ǷǺ㶨״̬��
!////��ʼϵ��
SUBROUTINE SideOutFlowInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    
    !TYPE(Com_element)::Com_element_main
    !TYPE(Condition)::Condition_0
    REAL::Waterdepth_flag
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation�洢ϵͳ��š����ʹ��롢�ֲ����
    !/////�ӳ�������������֮������ݽ���
    INTEGER::j_NUM  !Ԫ��ϵͳ���    
    INTEGER::Serial !��Ӧ��ˮ�ڵľֲ����
    REAL::h2,Q2,Q1,h1   !�ֱ��ʾѭ�������н��ڶ���ͳ��ڶ����ˮ�������
    INTEGER::nn
    INTEGER::i,k
    REAL::vvol  !Ԫ������
    
    !/////�ӳ����ڲ������漰���Ĳ���
    REAL::BB1,BB2   !��ˮ������
    REAL::Ar1,Ar2   !��ˮ���
    REAL::u1,u2 !��������
    REAL::aa,ee,dhh_S
    
    Ele_Rel=Com_element_main.Ele_Rel
    
    !///////��ˮ�ڷֶ���
    
    Com_element_main.n(j_NUM)=1!///��Զֻ��һ��
    !///////����ˮ�����
    h1=h2
    Do k=1,200,1
        aa=h1
        BB1=Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).BottomWidth+2*Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).SideSlopes*h1
        BB2=Com_element_main.Ca_cl(Ele_Rel(j_NUM+1,3)).BottomWidth+2*Com_element_main.Ca_cl(Ele_Rel(j_NUM+1,3)).SideSlopes*h2
        Ar1=(BB1+Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).BottomWidth)*h1/2.0
        Ar2=(BB2+Com_element_main.Ca_cl(Ele_Rel(j_NUM+1,3)).BottomWidth)*h2/2.0
        u1=Q1/Ar1
        u2=Q2/Ar2
        h1=h2+(u2**2-u1**2)/2.0/9.81 !����ˮͷ��
        IF(ABS(h1-aa)<2.22*2.718281828459**(-16))THEN   !����������ʵ�����������£������ж�
            EXIT
        END IF  
    END DO

    Waterdepth_flag=Waterdepth_flag+1
    Condition_0.h(Waterdepth_flag)=h1
    Condition_0.Q(Waterdepth_flag)=Q2
    vvol=(h1+h2)/2.0*((Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).BottomWidth+Com_element_main.Ca_cl(Ele_Rel(j_NUM+1,3)).BottomWidth)/2.0)*10.0 &         !�̶���ˮ�ڳ���Ϊ10m
         +Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).SideSlopes*h1*h2*10.0 &
         +1.0/3.0*(Com_element_main.Ca_cl(Ele_Rel(j_NUM-1,3)).SideSlopes*(h1-h2)**2)*10.0
END

!////////////////�Ǻ㶨��ϵ��
!////////////////�Ǻ㶨��ϵ��
Subroutine SideOutFlowDynamicCoeff(j_NUM,Serial,Qout,S2,T2,P2,V2,ele_vvol,Condition_cal)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_cal
    TYPE(Condition)::Condition_cal
    INTEGER::j_NUM  !Ԫ��ϵͳ���    
    INTEGER::Serial !��Ӧ�����ľֲ����
    REAL::Qout  !��ˮ����
    REAL::S2,T2,P2,V2   !׷��ϵ��
    REAL::C,D,E,F,G,M   !���Է������ϵ��
    REAL::Y1,Y2,Y3,Y4   !�м�ϵ��
    REAL::ele_vvol
    
    ele_vvol=0

    !/////��Ҫ�õ���ˮ��Ҫ��Ԥ����
    !//////���Է�����ϵ����ֵ
    C=0.0
    D=-Qout
    E=0.0
    G=0.0
    F=1.0
    M=0.0
    
    IF(Boundary_cal.Upsign==0)THEN
        !//////////�м�ϵ��
        Y1=Condition_cal.VV(j_NUM)+C
        Y2=F+E*Condition_cal.PP(j_NUM)
        Y3=D+Condition_cal.PP(j_NUM)
        Y4=M-E*Condition_cal.PP(j_NUM)

        !////////׷��ϵ��
        S2=(G*Y3-Y4)/(Y1*G+Y2)
        T2=(C*G-F)/(Y1*G+Y2)
        P2=Y3-Y1*S2
        V2=C-Y1*T2
    ELSE
        !////////׷��ϵ��
        S2=Qout
        T2=-1.0
        P2=Condition_cal.PP(j_NUM)+Condition_cal.VV(j_NUM)*(-Qout)
        V2=Condition_cal.VV(j_NUM)
    END IF
END 
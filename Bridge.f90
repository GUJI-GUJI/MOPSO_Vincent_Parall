!/////���� Ԫ������11
!////////�ֳ��������֣�һ���Ǻ㶨״̬�ģ�һ���ǷǺ㶨״̬��
!////��ʼϵ��
SUBROUTINE BridgeInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    
    !TYPE(Com_element)::Com_element_main
    !TYPE(Condition)::Condition_0
    REAL::Waterdepth_flag
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation�洢ϵͳ��š����ʹ��롢�ֲ����
    
    !/////�ӳ�������������֮������ݽ���
    INTEGER::j_NUM  !Ԫ��ϵͳ���    
    INTEGER::Serial !��Ӧ�������ľֲ����
    REAL::h2,Q2,Q1,h1   !�ֱ��ʾѭ�������н��ڶ���ͳ��ڶ����ˮ�������
    INTEGER::nn
    !/////�ӳ����ڲ������漰���Ĳ���
    REAL::BB1,BB2   !��ˮ������
    REAL::Ar_Brid   !�Ŷմ������
    REAL::Ar_NoBrid !�Ŷպ�����
    REAL::u_Brid    !�Ŷպ������
    REAL::alpha     !A'/A �Ŷ���ˮ��
    REAL::k         !�Ŷ���״ϵ��
    REAL::vvol
    
    Ele_Rel=Com_element_main.Ele_Rel
    
    !///////�����ֶ���

    Com_element_main.n(j_NUM)=1
    k=0.9   !///�Ŷ���״ϵ��

    !///////����ˮ�����
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


!////////////////�Ǻ㶨��ϵ��
Subroutine BridgeDynamicCoeff(j_NUM,Serial,Qo2,Qo1,Zo2,Zo1,S2,T2,P2,V2,ele_vvol,Condition_cal)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_cal
    TYPE(Condition)::Condition_cal

    INTEGER::j_NUM  !Ԫ��ϵͳ���    
    INTEGER::Serial !��Ӧ�������ľֲ����
    REAL::Qo2,Qo1,Zo2,Zo1 !�ֱ��ʾnʱ��j+1����������j����������j+1����ˮλ��j����ˮλ
    REAL::S2,T2,P2,V2   !׷��ϵ��
    REAL::C,D,E,F,G,M   !���Է������ϵ��
    REAL::Y1,Y2,Y3,Y4   !�м�ϵ��
    REAL::hh1,hh2   !ˮ��
    REAL::BB1,BB2   !ˮ����
    REAL::Ar1,Ar2,Ar_Brid      !��ˮ���
    REAL::L_wet     !ʪ��
    REAL::R !ˮ���뾶
    REAL::k !�Ŷ���״ϵ��
    REAL::alpha     !A'/A �Ŷ���ˮ��
    REAL::u1,u2 !��������
    REAL::ele_vvol  !Ԫ������

    !/////��Ҫ�õ���ˮ��Ҫ��Ԥ����
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

    !//////���Է�����ϵ����ֵ
    C=0.0
    D=0.0
    E=0.0
    G=0.0
    F=1.0
    M=-2*k*(k+10*((u2**2.0)/(2.0*9.81*hh2))-0.6)*(alpha+15.0*alpha**4.0)*((u2**2.0)/(2.0*9.81))

    
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
        P2=Condition_cal.PP(j_NUM)+M
        V2=Condition_cal.VV(j_NUM)
    ENDIF
END 
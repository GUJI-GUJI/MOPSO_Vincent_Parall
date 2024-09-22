!/////������ Ԫ������ 10
!////////�ֳ��������֣�һ���Ǻ㶨״̬�ģ�һ���ǷǺ㶨״̬��
!////��ʼϵ��
SUBROUTINE InvertedsiphonInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
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
    INTEGER::i
    INTEGER::nn
    REAL::vvol
    !/////�ӳ����ڲ������漰���Ĳ���
    REAL::BB1,BB2   !��ˮ������
    REAL::Ar_Insi   !�����ܵ����
    REAL::u_Insi    !�����ܵ�����
    REAL::L_wet !������ʪ��
    REAL::Hydr_R    !ˮ���뾶 !
    REAL::K_Insi    !����ģ��
    
    Ele_Rel=Com_element_main.Ele_Rel
    
    !///////�������ֶ���  
    Com_element_main.n(j_NUM)=1!///��Զֻ��һ��

    !///////����ˮ�����
    h1=h2
    Ar_Insi=Com_element_main.Insi_cl(Serial).ParallelNum*Com_element_main.Insi_cl(Serial).BottomWidth*Com_element_main.Insi_cl(Serial).Height
    BB1=Com_element_main.Insi_cl(Serial).ParallelNum*Com_element_main.Insi_cl(Serial).BottomWidth
    BB2=Com_element_main.Insi_cl(Serial).ParallelNum*Com_element_main.Insi_cl(Serial).BottomWidth
    L_wet=Com_element_main.Insi_cl(Serial).ParallelNum*(Com_element_main.Insi_cl(Serial).BottomWidth*2+2*Com_element_main.Insi_cl(Serial).Height)
    Hydr_R=Ar_Insi/L_wet
    u_Insi=Q2/Ar_Insi
    K_Insi=1.0/(Com_element_main.Insi_cl(Serial).Roughness)*Hydr_R**(1.0/6.0)*Ar_Insi*Hydr_R**(1.0/6.0)
    
    !���������غ㷽�̼���ˮ������ˮλ
    h1=Com_element_main.Insi_cl(Serial).OutInvertElev+h2-Com_element_main.Insi_cl(Serial).InInvertElev&
        &+Com_element_main.Insi_cl(Serial).Insi_coeff*Com_element_main.Insi_cl(Serial).Length/4.0/Hydr_R*u_Insi**2.0/2.0/9.81

    Waterdepth_flag=Waterdepth_flag+1
    Condition_0.h(Waterdepth_flag)=h1
    Condition_0.Q(Waterdepth_flag)=Q2
    vvol=(Com_element_main.Insi_cl(Serial).BottomWidth*Com_element_main.Insi_cl(Serial).ParallelNum-(Com_element_main.Insi_cl(Serial).ParallelNum-1)*&
        Com_element_main.Insi_cl(Serial).BottomWidth)*Com_element_main.Insi_cl(Serial).Height*Com_element_main.Insi_cl(Serial).Length
END

!////////////////�Ǻ㶨��ϵ��
!////////////////�Ǻ㶨��ϵ��
Subroutine InvertedsiphonDynamicCoeff(j_NUM,Serial,Qo2,Qo1,Zo2,Zo1,S2,T2,P2,V2,ele_vvol,Condition_cal)
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
    REAL::Ar1,Ar2,Ar_Insi      !��ˮ���
    REAL::L_wet     !ʪ��
    REAL::R !ˮ���뾶
    REAL::K !����ģ��
    REAL::C1,C2 !л��ϵ��
    REAL::u1,u2,u_Insi !��������
    REAL::ele_vvol

    !/////��Ҫ�õ���ˮ��Ҫ��Ԥ����
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
    !//////���Է�����ϵ����ֵ
    C=0.0
    D=0.0
    E=0.0
    G=0.0
    F=1.0

    !M=-(Insi_cl(Serial).Insi_coeff*u_Insi**2.0/2.0/GRAV+Qo1**2.0/K**2.0*Insi_cl(Serial).Length)
    M=-(Com_element_cal.Insi_cl(Serial).Insi_coeff*Com_element_cal.Insi_cl(Serial).Length/4.0/R*u_Insi**2.0/2.0/9.81)
    !M=-(Insi_cl(Serial).Insi_coeff)
    !M�ı��ʽ��Ҳ���ǵ�������ˮ����ʧ�����Ƕ��٣���ô��ʾ��ϵ���Ƕ��٣���һ��ֵ���о������⡣
    
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
    !/////����բ Ԫ������ 6
    !////////�ֳ��������֣�һ���Ǻ㶨״̬�ģ�һ���ǷǺ㶨״̬��
    !////��ʼϵ��
    SUBROUTINE ControlGateInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
    USE GLOBAL_VARIABLE

    IMPLICIT NONE
    !TYPE(Com_element)::Com_element_main
    !TYPE(Condition)::Condition_0
    REAL::Waterdepth_flag
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation�洢ϵͳ��š����ʹ��롢�ֲ����
    !/////�ӳ�������������֮������ݽ���
    INTEGER::j_NUM  !Ԫ��ϵͳ���
    INTEGER::Serial !��Ӧ����բ�ľֲ����
    REAL::h2,Q2,Q1,h1   !�ֱ��ʾѭ�������н��ڶ���ͳ��ڶ����ˮ�������
    INTEGER::i
    INTEGER::nn
    Real::Forebay_h !ǰ��ˮλ�Ϳ���
    REAL::vvol  !Ԫ������

    Forebay_h=Com_element_main.Cogate_cl(Serial).Initial_gate_Zup-Com_element_main.Cogate_cl(Serial).InvertElev
    !///////����բ�ֶ���
    Com_element_main.n(j_NUM)=1!///��Զֻ��һ��
    h1=Forebay_h
    Waterdepth_flag=Waterdepth_flag+1
    Condition_0.h(Waterdepth_flag)=h1
    Condition_0.Q(Waterdepth_flag)=Q2
    vvol=h2*Com_element_main.Cogate_cl(Serial).Lenth*Com_element_main.Cogate_cl(Serial).RUNNum*Com_element_main.Cogate_cl(Serial).SingleWidth
    END


    !////////////////�Ǻ㶨��ϵ��
    Subroutine ControlGateDynamicCoeff(j_NUM,Serial,Opde,Qo2,Qo1,Zo2,Zo1,S2,T2,P2,V2,ele_vvol,Condition_cal)
    USE GLOBAL_VARIABLE

    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_cal
    TYPE(Condition)::Condition_cal
    INTEGER::j_NUM  !Ԫ��ϵͳ���
    INTEGER::Serial !��Ӧ�����ľֲ����
    REAL::Qo2,Qo1,Zo2,Zo1 !�ֱ��ʾnʱ��j+1����������j����������j+1����ˮλ��j����ˮλ
    REAL::S2,T2,P2,V2   !׷��ϵ��
    REAL::C,D,E,F,G,M   !���Է������ϵ��
    REAL::Y1,Y2,Y3,Y4   !�м�ϵ��
    REAL::hh1,hh2   !ˮ��
    REAL::BB  !��բˮ�����
    REAL::Opde  !բ�ſ���
    REAL::CC    !�ۺ�ϵ��
    REAL::ele_vvol


    !/////��Ҫ�õ���ˮ��Ҫ��Ԥ����
    hh1=Zo1-Com_element_cal.Cogate_cl(Serial).InvertElev
    hh2=Zo2-Com_element_cal.Cogate_cl(Serial).InvertElev
    BB=Com_element_cal.Cogate_cl(Serial).RUNNum*Com_element_cal.Cogate_cl(Serial).SingleWidth
    ele_vvol=hh2*Com_element_cal.Cogate_cl(Serial).Lenth*Com_element_cal.Cogate_cl(Serial).RUNNum*Com_element_cal.Cogate_cl(Serial).SingleWidth
    if ((hh1<hh2)) then
        write(*,*)
        write(*,*)
        write(*,*)"=======WARNING��ControlGate.f90�����Ǻ㶨�����㣩====="
        write(*,*) "��",j_NUM,"������Ԫ����"
        write(*,*) "����բբǰˮλ����բ��ˮλ"
        write(*,*) "˵����һʱ�̼����ˮ��������Զ����Զ�����"
        hh2=hh1-0.05
    end if

    !!�̶�ϵ��
    CC=(Com_element_cal.Cogate_cl(Serial).Discoeff)*BB*Opde*SQRT(2.0*9.81)   !ģ���� �̶�ϵ��
    !!բ�ſ���-����ϵ��
    !CC=(Com_element_cal.Cogate_cl(Serial).G_A*Opde**2+Com_element_cal.Cogate_cl(Serial).G_B*Opde+Com_element_cal.Cogate_cl(Serial).G_C)*BB*Opde*SQRT(2.0*9.81)  !��ϵ��
    !!բ�ſ���/ˮͷ��-����ϵ��
    !CC=(Com_element_cal.Cogate_cl(Serial).G_A*(Opde/(Zo1-Zo2))**2+Com_element_cal.Cogate_cl(Serial).G_B*(Opde/(Zo1-Zo2))+Com_element_cal.Cogate_cl(Serial).G_C)*BB*Opde*SQRT(2.0*9.81)




    !IF(Boundary_cal.Upsign==0)THEN
    !    !//////////�м�ϵ��
    !    Y1=Condition_cal.VV(j_NUM)+C
    !    Y2=F+E*Condition_cal.VV(j_NUM)
    !    Y3=D+Condition_cal.PP(j_NUM)
    !    Y4=M-E*Condition_cal.PP(j_NUM)
    !
    !    !////////׷��ϵ��
    !    S2=(G*Y3-Y4)/(Y1*G+Y2)
    !    T2=(C*G-F)/(Y1*G+Y2)
    !    P2=Y3-Y1*S2
    !    V2=C-Y1*T2
    !ELSE
    !    !////////׷��ϵ��
    !    S2=0.0
    !    T2=-1.0
    !    P2=M/F+Condition_cal.PP(j_NUM)
    !    V2=Condition_cal.VV(j_NUM)+1.0/F
    !ENDIF


    !׷��ϵ��
    IF (Opde<=MAX(hh1,hh2))Then !բ�Ž���
        !����բǰˮλ����С
        IF(abs(hh1-hh2)<0.01)Then
            hh1=hh1+0.02
            hh2=hh2-0.02
        ENDIF

        !����բǰˮλ����բ��ˮλ
        IF(hh1<hh2)THEN
            hh1=hh2+0.05
            Write(*,'(*(G0.5))') "����բ������բǰˮλ����բ��ˮλ���������ǿ����������ע����֤����ɿ��ԣ�"
        ENDIF

        !//////���Է�����ϵ����ֵ(����������û������
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


    ELSEIF(Opde>MAX(hh1,hh2))Then!բ��ȫ��
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
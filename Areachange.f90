!/////����� Ԫ������7
!////////�ֳ��������֣�һ���Ǻ㶨״̬�ģ�һ���ǷǺ㶨״̬��
!////��ʼϵ��
SUBROUTINE AreachangeInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Config)::Config_main
    !!TYPE(Com_element)::Com_element_main
    !TYPE(Condition)::Condition_0
    REAL::GRAV
    REAL::Waterdepth_flag
    
    !/////�ӳ�������������֮������ݽ���
    INTEGER::j_NUM  !Ԫ��ϵͳ���
    INTEGER::Serial !��Ӧ����εľֲ����
    REAL::h2,Q2,Q1,h1   !�ֱ��ʾѭ�������г��ڶ���ͽ��ڶ����ˮ�������
    INTEGER::nn
    INTEGER::i,k
    REAL::vvol   !Ԫ������
    
    !/////�ӳ����ڲ������漰���Ĳ���
    REAL::BB1,BB2   !��ˮ������
    REAL::Ar1,Ar2   !��ˮ���
    REAL::u1,u2 !��������
    REAL::aa,ee,dhh_S
    REAL::delta_elev   !�����ڸ̲߳�
    REAL::BW_in    !������ڵ׿�
    REAL::BW_out    !������ڵ׿�
    
    !!!!!!////�ر�ע�⣺������Ҫ��������ۣ�Ҫ�������ǰ��ӵ���ʲô���п�����������������Ҳ�п����������ӵ�����
    
    GRAV=Config_main.Coeff_set.GRAV
    
    !///////����ηֶ���

    Com_element_main.n(j_NUM)=1
    !//////���γ�ʼˮ��
    IF(Com_element_main.Ele_Rel(j_NUM-1,2)==5.AND.Com_element_main.Ele_Rel(j_NUM+1,2)==5)THEN
        h1=h2
        Do k=1,200,1
            aa=h1
            delta_elev=Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).OutInvertElev-Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).InInvertElev
            BW_in=Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).BottomWidth
            BW_out=Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).BottomWidth
            BB1=Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).BottomWidth+2*Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).SideSlopes*h1
            BB2=Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).BottomWidth+2*Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).SideSlopes*h2
            Ar1=(BB1+Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).BottomWidth)*h1/2.0
            Ar2=(BB2+Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).BottomWidth)*h2/2.0
            u1=Q1/Ar1
            u2=Q2/Ar2
            h1=h2+Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).InInvertElev-Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).OutInvertElev& !�߳�
                &+(u2**2-u1**2)/2.0/GRAV& !����ˮͷ��
                &+Com_element_main.Arch_cl(Serial).Locallosscoeff*(u2/2.0+u1/2.0)**2.0/2.0/GRAV&   !�ֲ���ʧ
                &+Com_element_main.Arch_cl(Serial).Additionlosscoeff*abs(u1-u2)**2.0/2.0/GRAV !�س���ʧ
            IF(ABS(h1-aa)<2.22*2.718281828459**(-16))THEN   !����������ʵ�����������£������ж�
                EXIT
            END IF  
        END DO
    ELSEIF(Com_element_main.Ele_Rel(j_NUM-1,2)==5.AND.Com_element_main.Ele_Rel(j_NUM+1,2)==10)THEN    !��������ϵ������Ҫ�úÿ���һ��
        h1=h2
        Do k=1,200,1
            aa=h1
            delta_elev=Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).OutInvertElev-Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).InInvertElev
            BW_in=Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).BottomWidth
            BW_out=Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).BottomWidth*Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).ParallelNum
            BB1=Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).BottomWidth+2.0*Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).SideSlopes*h1
            BB2=Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).ParallelNum*Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).BottomWidth
            Ar1=(BB1+Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).BottomWidth)*h1/2.0
            Ar2=BB2*h2
            u1=Q1/Ar1
            u2=Q2/Ar2
            h1=h2+Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).InInvertElev-Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).OutInvertElev& !�߳�
                &+(u2**2-u1**2)/2.0/GRAV& !����ˮͷ��
                &+Com_element_main.Arch_cl(Serial).Locallosscoeff*(u2/2.0+u1/2.0)**2.0/2.0/GRAV&   !�ֲ���ʧ
                &+Com_element_main.Arch_cl(Serial).Additionlosscoeff*abs(u1-u2)**2.0/2.0/GRAV !�س���ʧ
            IF(ABS(h1-aa)<2.22*2.718281828459**(-16))THEN   !����������ʵ�����������£������ж�
                EXIT
            END IF  
        END DO
    ELSEIF(Com_element_main.Ele_Rel(j_NUM-1,2)==10.AND.Com_element_main.Ele_Rel(j_NUM+1,2)==5)THEN    !��������ϵ������Ҫ�úÿ���һ��
        h1=h2
        Do k=1,200,1
            aa=h1
            delta_elev=Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).OutInvertElev-Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).InInvertElev
            BW_in=Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).BottomWidth*Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).ParallelNum
            BW_out=Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).BottomWidth
            BB1=Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).ParallelNum*Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).BottomWidth
            BB2=Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).BottomWidth+2*Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).SideSlopes*h2
            Ar1=BB1*h1
            Ar2=(BB2+Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).BottomWidth)*h2/2.0
            u1=Q1/Ar1
            u2=Q2/Ar2
            h1=h2+Com_element_main.Ca_cl(Com_element_main.Ele_Rel(j_NUM+1,3)).InInvertElev-Com_element_main.Insi_cl(Com_element_main.Ele_Rel(j_NUM-1,3)).OutInvertElev& !�߳�
               &+(u2**2-u1**2)/2.0/GRAV& !����ˮͷ��
               &+Com_element_main.Arch_cl(Serial).Locallosscoeff*(u2/2.0+u1/2.0)**2.0/2.0/GRAV&   !�ֲ���ʧ
               &+Com_element_main.Arch_cl(Serial).Additionlosscoeff*abs(u1-u2)**2.0/2.0/GRAV !�س���ʧ
            IF(ABS(h1-aa)<2.22*2.718281828459**(-16))THEN   !����������ʵ�����������£������ж�
                EXIT
            END IF  
        END DO
    END IF
    Waterdepth_flag=Waterdepth_flag+1
    Condition_0.h(Waterdepth_flag)=h1
    Condition_0.Q(Waterdepth_flag)=Q2
    vvol=1.0/6.0*BW_in*(2*h2+2*delta_elev+h2)*1000.0*(Com_element_main.Arch_cl(Serial).Outcoord-Com_element_main.Arch_cl(Serial).Incoord) &
            &+1.0/6.0*BW_out*(2*h2+h2+delta_elev)*1000.0*(Com_element_main.Arch_cl(Serial).Outcoord-Com_element_main.Arch_cl(Serial).Incoord)&
            &+Com_element_main.Arch_cl(Serial).SideSlopes*(h2+delta_elev)*h2*1000.0*(Com_element_main.Arch_cl(Serial).Outcoord-Com_element_main.Arch_cl(Serial).Incoord) &
            &+1.0/3.0*Com_element_main.Arch_cl(Serial).SideSlopes*(delta_elev**2)*1000.0*(Com_element_main.Arch_cl(Serial).Outcoord-Com_element_main.Arch_cl(Serial).Incoord)
END

!///////�Ǻ㶨������
Subroutine AreachangeDynamicCoeff(j_NUM,Serial,Qo2,Qo1,Zo2,Zo1,S2,T2,P2,V2,ele_vvol,Condition_cal)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_cal
    !TYPE(Config)::Config_main
    TYPE(Condition)::Condition_cal
    REAL::GRAV
    
    
    INTEGER::j_NUM  !Ԫ��ϵͳ���    
    INTEGER::Serial !��Ӧ�����ľֲ����
    REAL::Qo2,Qo1,Zo2,Zo1 !�ֱ��ʾnʱ��j+1����������j����������j+1����ˮλ��j����ˮλ
    REAL::S2,T2,P2,V2   !׷��ϵ��
    REAL::C,D,E,F,G,M   !���Է������ϵ��
    REAL::Y1,Y2,Y3,Y4   !�м�ϵ��
    REAL::hh1,hh2   !ˮ��
    REAL::BB1,BB2   !���
    REAL::delta_elev  !�����ڸ̲߳�
    REAL::BW_in    !������ڵ׿�
    REAL::BW_out    !������ڵ׿�
    REAL::Ar1,Ar2      !��ˮ���
    REAL::u1,u2 !��������
    REAL::ele_vvol
    

    GRAV=Config_main.Coeff_set.GRAV

    !/////��Ҫ�õ���ˮ��Ҫ��Ԥ����
    IF(Com_element_cal.Ele_Rel(j_NUM-1,2)==5.AND.Com_element_cal.Ele_Rel(j_NUM+1,2)==5)THEN
        hh1=Zo1-Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).OutInvertElev
        hh2=Zo2-Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).InInvertElev
        delta_elev=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).OutInvertElev-Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).InInvertElev
        BW_in=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).BottomWidth
        BW_out=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).BottomWidth
        BB1=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).BottomWidth+2*Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).SideSlopes*hh1
        BB2=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).BottomWidth+2*Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).SideSlopes*hh2
        Ar1=(BB1+Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).BottomWidth)*hh1/2.0
        Ar2=(BB2+Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).BottomWidth)*hh2/2.0
    ELSEIF(Com_element_cal.Ele_Rel(j_NUM-1,2)==5.AND.Com_element_cal.Ele_Rel(j_NUM+1,2)==10)THEN
        hh1=Zo1-Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).OutInvertElev
        hh2=Zo2-Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).InInvertElev
        delta_elev=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).OutInvertElev-Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).InInvertElev
        BW_in=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).BottomWidth
        BW_out=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).BottomWidth*Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).ParallelNum
        BB1=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).BottomWidth+2.0*Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).SideSlopes*hh1
        BB2=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).ParallelNum*Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).BottomWidth
        Ar1=(BB1+Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).BottomWidth)*hh1/2.0
        Ar2=BB2*hh2
    ELSEIF(Com_element_cal.Ele_Rel(j_NUM-1,2)==10.AND.Com_element_cal.Ele_Rel(j_NUM+1,2)==5)THEN
        hh1=Zo1-Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).OutInvertElev
        hh2=Zo2-Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).InInvertElev
        delta_elev=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).OutInvertElev-Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).InInvertElev
        BW_in=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).BottomWidth*Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).ParallelNum
        BW_out=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).BottomWidth
        BB1=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).ParallelNum*Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j_NUM-1,3)).BottomWidth
        BB2=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).BottomWidth+2*Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).SideSlopes*hh2
        Ar1=BB1*hh1
        Ar2=(BB2+Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j_NUM+1,3)).BottomWidth)*hh2/2.0
    END IF
    
    ele_vvol=1.0/6.0*BW_in*(2*hh2+2*delta_elev+hh2)*1000.0*(Com_element_cal.Arch_cl(Serial).Outcoord-Com_element_cal.Arch_cl(Serial).Incoord) &
            &+1.0/6.0*BW_out*(2*hh2+hh2+delta_elev)*1000.0*(Com_element_cal.Arch_cl(Serial).Outcoord-Com_element_cal.Arch_cl(Serial).Incoord)&
            &+Com_element_cal.Arch_cl(Serial).SideSlopes*(hh2+delta_elev)*hh2*1000.0*(Com_element_cal.Arch_cl(Serial).Outcoord-Com_element_cal.Arch_cl(Serial).Incoord) &
            &+1.0/3.0*(Com_element_cal.Arch_cl(Serial).SideSlopes*delta_elev**2)*1000.0*(Com_element_cal.Arch_cl(Serial).Outcoord-Com_element_cal.Arch_cl(Serial).Incoord)
    
    !//////�����������
    u1=Qo1/Ar1
    u2=Qo2/Ar2

    !//////���Է�����ϵ����ֵ
    C=0.0
    D=0.0
    E=0.0
    G=0.0
    F=1.0
    M=-((u2**2-u1**2)/2.0/GRAV+Com_element_cal.Arch_cl(Serial).Locallosscoeff*(u2/2.0+u1/2.0)**2.0/2.0/GRAV+Com_element_cal.Arch_cl(Serial).Additionlosscoeff*abs(u1-u2)**2.0/2.0/GRAV)
    
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
    END IF
END 
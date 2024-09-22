!////������Ԫ������ 5
!////////�ֳ��������֣�һ���Ǻ㶨״̬�ģ�һ���ǷǺ㶨״̬��
!////��ʼϵ��
SUBROUTINE CanalInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_main
    !TYPE(Condition)::Condition_0
    !TYPE(Config)::Config_main
    REAL::Psi
    REAL::PsiR
    REAL::GRAV
    REAL::Waterdepth_flag
    
    !/////�ӳ�������������֮������ݽ���
    INTEGER::j_NUM  !Ԫ��ϵͳ���
    INTEGER::Serial !��Ӧ�����ľֲ����
    REAL::h2,Q2,Q1,h1   !�ֱ��ʾѭ�������н��ڶ���ͳ��ڶ����ˮ�������
    INTEGER::nn  !Ԫ���ֶθ���
    INTEGER::j,i
    REAL::aa,ee,dhh_S
    REAL::vvol  !Ԫ������
    REAL::DeltaX    !�ռ䲽��

    !//////������ˮ�����
    REAL,ALLOCATABLE::QQ_S(:)   !Ԫ������������
    REAL,ALLOCATABLE::hh_S(:)   !Ԫ��������ˮλ

    !/////�ӳ����ڲ������漰���Ĳ���
    REAL::Len   !Ԫ������ !һ�㲻�ùؼ��������� wjb0415
    REAL::S0    !Ԫ������

    !//////////�洢����СԪ���Ĳ�����n��Ԫ��
    REAL,ALLOCATABLE::BottomWidth_s(:)    !��СԪ���׿�
    REAL,ALLOCATABLE::SideSlopes_s(:) !��СԪ������ϵ��
    REAL,ALLOCATABLE::Incoord_s(:) !��СԪ��������
    REAL,ALLOCATABLE::Outcoord_s(:) !��СԪ���������
    REAL,ALLOCATABLE::InInvertElev_s(:)   !��СԪ����ڵ׸߳�
    REAL,ALLOCATABLE::OutInvertElev_s(:)   !��СԪ�����ڵ׸߳�
    REAL,ALLOCATABLE::InTopElev_s(:)   !��СԪ����ڶ��߳�
    REAL,ALLOCATABLE::OutTopElev_s(:)   !��СԪ�����ڶ��߳�
    
    !/////////�洢����������ݣ�n+1������
    REAL,ALLOCATABLE::BB(:) !������ˮ����
    REAL,ALLOCATABLE::Ar(:) !�������ˮ���
    REAL,ALLOCATABLE::Sf(:) !������Ħ���¶�
    REAL,ALLOCATABLE::KK(:) !����������ģ��
    REAL,ALLOCATABLE::L_wet(:)  !������ʪ�ܳ���
    REAL,ALLOCATABLE::RR(:) !������ˮ���뾶
    !------------------
    INTEGER :: CroSec_TYPE
    !-----------------------------

    Psi=Config_main.Coeff_set.Psi
    PsiR=Config_main.Coeff_set.PsiR
    GRAV=Config_main.Coeff_set.GRAV
    
    DeltaX=Boundary_cal.Step*SQRT(GRAV*h2)
    !WRITE(*,*)DeltaX
    Len=(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord)*1000  

    Com_element_main.n(j_NUM)=int(Len/DeltaX)
    IF(Com_element_main.n(j_NUM)<1)THEN
        Com_element_main.n(j_NUM)=1
    END IF

    if (Com_element_main.Ca_cl(Serial).CroSecTYPE==2) Com_element_main.n(j_NUM)=1!#TU�ӣ����������Ļ����ֶμ���
    nn=Com_element_main.n(j_NUM)
    DeltaX=Len/nn    !�ռ䲽��
    ALLOCATE(BottomWidth_s(nn))  
    ALLOCATE(SideSlopes_s(nn))
    ALLOCATE(Incoord_s(nn))
    ALLOCATE(Outcoord_s(nn)) 
    ALLOCATE(InInvertElev_s(nn))   
    ALLOCATE(OutInvertElev_s(nn))  
    ALLOCATE(InTopElev_s(nn))   
    ALLOCATE(OutTopElev_s(nn)) 
    ALLOCATE(BB(nn+1))
    ALLOCATE(Ar(nn+1))
    ALLOCATE(Sf(nn+1))
    ALLOCATE(KK(nn+1))
    ALLOCATE(QQ_S(nn+1))
    ALLOCATE(hh_S(nn+1))
    ALLOCATE(L_wet(nn+1))
    ALLOCATE(RR(nn+1))

    DO j=1,nn
        BottomWidth_s(j)=Com_element_main.Ca_cl(Serial).BottomWidth
        SideSlopes_s(j)=Com_element_main.Ca_cl(Serial).SideSlopes
        Incoord_s(j)=Com_element_main.Ca_cl(Serial).Incoord+(j-1)*(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord)/nn
        Outcoord_s(j)=Com_element_main.Ca_cl(Serial).Incoord+j*(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord)/nn
        InInvertElev_s(j)=Com_element_main.Ca_cl(Serial).InInvertElev+(j-1)*(Com_element_main.Ca_cl(Serial).OutInvertElev-Com_element_main.Ca_cl(Serial).InInvertElev)/nn
        OutInvertElev_s(j)=Com_element_main.Ca_cl(Serial).InInvertElev+j*(Com_element_main.Ca_cl(Serial).OutInvertElev-Com_element_main.Ca_cl(Serial).InInvertElev)/nn
        InTopElev_s(j)=Com_element_main.Ca_cl(Serial).InTopElev+(j-1)*(Com_element_main.Ca_cl(Serial).OutTopElev-Com_element_main.Ca_cl(Serial).InTopElev)/nn
        OutTopElev_s(j)=Com_element_main.Ca_cl(Serial).InTopElev+j*(Com_element_main.Ca_cl(Serial).OutTopElev-Com_element_main.Ca_cl(Serial).InTopElev)/nn
    END DO
    
    S0=(Com_element_main.Ca_cl(Serial).InInvertElev-Com_element_main.Ca_cl(Serial).OutInvertElev)/((Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord)*1000)
    DO j=1,nn+1
        QQ_S(j)=Q1    
    END DO
    !--------------------wjb0414
    CroSec_TYPE=Com_element_main.Ca_cl(Serial).CroSecTYPE
    !-------------------------
    hh_S(nn+1)=h2
    hh_S(nn)=hh_S(nn+1)

   !������㣬���ֶ������ͣ�wjb0414 ����һ���������Ͳ����Ͷ�������Ӻ���
    DO j=nn,1,-1
        hh_S(j)=hh_S(j+1)
        DO i=1,200,1 !////200�����ѭ����
            !j+1����
            CALL Abr(Serial,CroSec_TYPE,hh_S(j+1),BottomWidth_s(j),SideSlopes_s(j),OutInvertElev_s(j),BB(j+1),Ar(j+1),L_wet(j+1),1,0)  !#TU
            !ˮ���뾶
            RR(j+1)=Ar(j+1)/L_wet(j+1)
            !����ģ��K=CAsqrt(R)
            KK(j+1)=1/Com_element_main.Ca_cl(Serial).Roughness*Ar(j+1)**(5.0/3.0)*L_wet(j+1)**(-2.0/3.0) !ͬ��ʱ�Ѳ������ó���ͬ�ľͿ����ˡ�
            Sf(j+1)=QQ_S(j+1)*ABS(QQ_S(j+1))/KK(j+1)**2.0
            !j����
            CALL Abr(Serial,CroSec_TYPE,hh_S(j),BottomWidth_s(j),SideSlopes_s(j),inInvertElev_s(j),BB(j),Ar(j),L_wet(j),0,0)!#TU

            RR(j)=Ar(j)/L_wet(j)
            KK(j)=1/Com_element_main.Ca_cl(Serial).Roughness*Ar(j)**(5.0/3.0)*L_wet(j)**(-2.0/3.0)!ͬ��ʱ�Ѳ������ó���ͬ�ľͿ����ˡ�
            Sf(j)=QQ_S(j)*ABS(QQ_S(j))/KK(j)**2.0

            !��������
            !Tu BottomWidth_s(j)+2*SideSlopes_s(j)*hh_S(j)��һ������ΪBB(j)�����ⲻ�������׿�Ϊ�̶߳Ե����⣻����ϵ��������Ӱ�첻���Ȳ���
            aa=Psi/DeltaX*QQ_S(j)**2*BB(j)/(Ar(j)**3)-GRAV*Psi/DeltaX-2.0*Psi*(1.0-PsiR)*GRAV*Sf(j)/KK(j)*KK(j)/Ar(j)*(5.0/3.0*BB(j)-2.0/3.0*RR(j)*(2.0*SQRT(1.0+SideSlopes_s(j)**2)))!#TU��
            ee=-Psi/DeltaX/2.0*((QQ_S(j+1)/Ar(j+1))**2.0-(QQ_S(j)/Ar(j))**2.0)-(1.0-Psi)/DeltaX/2.0*((QQ_S(j+1)/Ar(j+1))**2.0-(QQ_S(j)/Ar(j))**2.0)-Psi*GRAV/DeltaX*(hh_S(j+1)-hh_S(j))-(1.0-Psi)*GRAV/DeltaX*(hh_S(j+1)-hh_S(j))-Psi*GRAV*(PsiR*Sf(j+1)+(1.0-PsiR)*Sf(j))-(1.0-Psi)*GRAV*(PsiR*Sf(j+1)+(1.0-PsiR)*Sf(j))+GRAV*S0

            !��������
            dhh_S=ee/aa
            hh_S(j)=hh_S(j)+dhh_S
            !TU���������һ����һ��aa���ӽ�0��С��0��
            !��ee����0����ôhh_S�Ǹ��ģ�����һϵ�ж�����

            !#TU��
            if (hh_S(j)<0) then
                write(*,*)"==============warning���㶨��������======================"
                write(*,*) "�ڵ�",j_NUM,"����Ԫ����������"
                write(*,*) "��",nn,"��С������"
                write(*,*) "�㶨��ˮ�����������������ˮλС��0�������޸�Ϊ3��������"
                write(*,*) "����˾�����ֺܶ�������������ȥ�����Գ����޸�3Ϊ����ֵ"
                write(*,*)"======================================================="
                hh_S(j)=3
            end if

            IF(ABS(dhh_S)<2.22*2.718281828459**(-16))THEN
                EXIT
            END IF
        END DO
        Waterdepth_flag=Waterdepth_flag+1
        Condition_0.h(Waterdepth_flag)=hh_S(j)
        Condition_0.Q(Waterdepth_flag)=QQ_S(j)
    END DO
    
    h1=hh_S(1)
    vvol=(h1+h2)/2.0*Com_element_main.Ca_cl(Serial).BottomWidth*1000.0*(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord) &
         &+Com_element_main.Ca_cl(Serial).SideSlopes*h1*h2*1000.0*(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord) &
         &+1.0/3.0*(Com_element_main.Ca_cl(Serial).SideSlopes*(h1-h2)**2)*1000.0*(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord)
    DEALLOCATE(BottomWidth_s)    
    DEALLOCATE(SideSlopes_s)
    DEALLOCATE(Incoord_s)
    DEALLOCATE(Outcoord_s) 
    DEALLOCATE(InInvertElev_s)   
    DEALLOCATE(OutInvertElev_s)  
    DEALLOCATE(InTopElev_s)   
    DEALLOCATE(OutTopElev_s)
    DEALLOCATE(BB)
    DEALLOCATE(Ar)
    DEALLOCATE(Sf)
    DEALLOCATE(KK)
    DEALLOCATE(QQ_S)
    DEALLOCATE(hh_S)
    DEALLOCATE(L_wet)
    DEALLOCATE(RR)
END

!////////////////�Ǻ㶨��ϵ��
Subroutine CanalDynamicCoeff(j_NUM,Serial,Qo2,Qo1,Zo2,Zo1,S2,T2,P2,V2,ele_vvol,Condition_cal)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_cal
    TYPE(Condition)::Condition_cal
    !TYPE(Config)::Config_main
    REAL::Theta
    REAL::GRAV
    REAL::RoughnessCoeff

    
    INTEGER::j_NUM  !Ԫ��ϵͳ���    
    INTEGER::Serial !��Ӧ�����ľֲ����
    REAL::Qo2,Qo1,Zo2,Zo1 !�ֱ��ʾnʱ��j+1����������j����������j+1����ˮλ��j����ˮλ
    REAL::S2,T2,P2,V2   !׷��ϵ��
    REAL::C,D,E,F,G,M   !���Է������ϵ��
    REAL::Y1,Y2,Y3,Y4   !�м�ϵ��
    REAL::hh1,hh2   !ˮ��
    REAL::BB1,BB2   !ˮ����
    REAL::Ar1,Ar2      !��ˮ���
    REAL::L_wet1,L_wet2     !ʪ��
    REAL::R1,R2 !ˮ���뾶
    REAL::C1,C2 !л��ϵ��
    REAL::u1,u2 !��������
    integer ::k             !#Tu����ʱ�Σ����ڱ���
    REAL::ele_vvol  !Ԫ������
    REAL::DeltaX    !�ռ䲽��
    
    Theta=Config_main.Coeff_set.Theta
    GRAV=Config_main.Coeff_set.GRAV
    RoughnessCoeff=Config_main.Coeff_set.RoughnessCoeff


    !/////��Ҫ�õ���ˮ��Ҫ��Ԥ����
    hh1=Zo1-Com_element_cal.Ca_cl(Serial).InInvertElev
    hh2=Zo2-Com_element_cal.Ca_cl(Serial).OutInvertElev
    !#Tu-----------------------------------------���¾����޸�
    if ((hh1<0).or.(hh2<0)) then
        write(*,*)
        write(*,*)
        write(*,*)"==============warning���Ǻ㶨�����㣩==============="
        write(*,*) "�ڵ�",k,"ʱ�̼�����"
        write(*,*) "��",j_NUM,"������Ԫ������������"
        write(*,*) "�Ǻ㶨������ˮ��С��0"
        write(*,*) "˵����һʱ�̼����ˮ�����"
        write(*,*)"=================================================="
    end if

    CALL Abr(Serial,Com_element_cal.Ca_cl(Serial).CroSecTYPE,hh2,Com_element_cal.Ca_cl(Serial).BottomWidth,Com_element_cal.Ca_cl(Serial).SideSlopes,Com_element_cal.Ca_cl(Serial).OutInvertElev,BB2,Ar2,L_wet2,1,1)!TU#@�׿����� ��������ζ��棬һ������������һ���ģ�����ǲ�����������ò���
    CALL Abr(Serial,Com_element_cal.Ca_cl(Serial).CroSecTYPE,hh1,Com_element_cal.Ca_cl(Serial).BottomWidth,Com_element_cal.Ca_cl(Serial).SideSlopes,Com_element_cal.Ca_cl(Serial).InInvertElev,BB1,Ar1,L_wet1,0,1)!TU#�׸߳����෴��ֻ�в���������õĵ��������������,��������ǿ��Ʒ�

    if ((Ar1<0.001).or.(Ar2<0.001)) then
        write(*,*)
        write(*,*)
        write(*,*)"==============warning���Ǻ㶨�����㣩==============="
        write(*,*) "�ڵ�",k,"ʱ�̼�����"
        write(*,*) "��",j_NUM,"������Ԫ������������"
        write(*,*) "���Ϊ0"
        write(*,*) "˵��������"
        write(*,*)"=================================================="
    end if

    !BB1=Com_element_cal.Ca_cl(Serial).BottomWidth+2.0*Com_element_cal.Ca_cl(Serial).SideSlopes*hh1
    !BB2=Com_element_cal.Ca_cl(Serial).BottomWidth+2.0*Com_element_cal.Ca_cl(Serial).SideSlopes*hh2
    DeltaX=(Com_element_cal.Ca_cl(Serial).Outcoord-Com_element_cal.Ca_cl(Serial).Incoord)*1000.0
    !Ar1=(BB1+Com_element_cal.Ca_cl(Serial).BottomWidth)*hh1/2.0
    !Ar2=(BB2+Com_element_cal.Ca_cl(Serial).BottomWidth)*hh2/2.0
    !L_wet1=Com_element_cal.Ca_cl(Serial).BottomWidth+2.0*hh1*SQRT(1+(Com_element_cal.Ca_cl(Serial).SideSlopes)**2.0)
    !L_wet2=Com_element_cal.Ca_cl(Serial).BottomWidth+2.0*hh2*SQRT(1+(Com_element_cal.Ca_cl(Serial).SideSlopes)**2.0)
    R1=Ar1/L_wet1
    R2=Ar2/L_wet2!#@TU
    u1=Qo1/Ar1
    u2=Qo2/Ar2
    !#Tu-----------------------------------------�޸ĵ��˽���
    
    ele_vvol=(hh1+hh2)/2.0*Com_element_cal.Ca_cl(Serial).BottomWidth*1000.0*(Com_element_cal.Ca_cl(Serial).Outcoord-Com_element_cal.Ca_cl(Serial).Incoord) &
         +Com_element_cal.Ca_cl(Serial).SideSlopes*hh1*hh2*1000.0*(Com_element_cal.Ca_cl(Serial).Outcoord-Com_element_cal.Ca_cl(Serial).Incoord) &
         +1.0/3.0*(Com_element_cal.Ca_cl(Serial).SideSlopes*(hh1-hh2)**2)*1000.0*(Com_element_cal.Ca_cl(Serial).Outcoord-Com_element_cal.Ca_cl(Serial).Incoord)
    
    IF(Config_main.Assimilate_set.Control_sign==0)THEN     !ģ����
        C1=1.0/Com_element_cal.Ca_cl(Serial).Roughness*R1**(1.0/6.0)
        C2=1.0/Com_element_cal.Ca_cl(Serial).Roughness*R2**(1.0/6.0)    
    ELSE                        !ͬ����
        C1=1.0/RoughnessCoeff*R1**(1.0/6.0)
        C2=1.0/RoughnessCoeff*R2**(1.0/6.0)
    END IF

    !//////���Է�����ϵ����ֵ
    C=1.0/2.0*(BB1+BB2)*DeltaX/2.0/Boundary_cal.Step/Theta
    D=-(1.0-Theta)/Theta*(Qo2-Qo1)+C*(Zo2+Zo1)
    E=DeltaX/2.0/Theta/Boundary_cal.Step-u1+GRAV*ABS(u1)/2.0/Theta/C1**2/R1*DeltaX
    G=DeltaX/2.0/Theta/Boundary_cal.Step+u2+GRAV*ABS(u2)/2.0/Theta/C2**2/R2*DeltaX
    F=GRAV*1.0/2.0*(Ar1+Ar2)
    M=DeltaX/2.0/Theta/Boundary_cal.Step*(Qo2+Qo1)-(1-Theta)/Theta*(u2*Qo2-u1*Qo1)&
      -(1-Theta)/Theta*GRAV*1.0/2.0*(Ar1+Ar2)*(Zo2-Zo1)
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
        !//////////�м�ϵ��
        Y1=D-C*Condition_cal.PP(j_NUM)
        Y2=M+F*Condition_cal.PP(j_NUM)
        Y3=1.0+C*Condition_cal.VV(j_NUM)
        Y4=E+F*Condition_cal.VV(j_NUM)

        !////////׷��ϵ��
        S2=(C*Y2-F*Y1)/(F*Y3+C*Y4)
        T2=(C*G-F)/(F*Y3+C*Y4)
        P2=(Y1+Y3*S2)/C
        V2=(Y3*T2+1.0)/C
    END IF
END 
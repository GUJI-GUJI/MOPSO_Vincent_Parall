    SUBROUTINE Dynamic_Cal(Boundary_temp,pollu_Lenth,arrive_time)  !�Ǻ㶨������
    USE GLOBAL_VARIABLE

    IMPLICIT NONE

    TYPE(Boundary)::Boundary_temp


    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_main,Com_element_cal
    !TYPE(Config)::Config_main

    TYPE(Dynamic_Result)::Dynamic_Result_cal
    !TYPE(Condition)::Condition_0,Condition_cal
    TYPE(Condition)::Condition_cal
    INTEGER::NUMMember
    INTEGER,ALLOCATABLE::DM_POSTION(:)  !����բ+��β����λ��
    INTEGER,ALLOCATABLE::SG_POSTION(:)  !����բλ��
    REAL::m_coeff    !ĩ��բ�Ź�բ����ϵ������
    REAL::PPQ, VVQ  !ˮλ~������ϵ����ϵ��
    INTEGER::Canal_NUM
    INTEGER::c_element    !��������ռ䲽����������Ԫ���ĸ���

    INTEGER::i,j,k,mm,kk,LLL
    real::c11,c22,c21,hh
    REAL::zzaverage,QQaverage
    INTEGER::si_n   !�Ż��洢ʱ��ѭ������
    REAL::random    !�����
    DOUBLE PRECISION r0
    Integer::pos_sg,DM_sg !��¼����բλ��
    INTEGER::ele_sg    !��¼��Ԫ����λ�ã������Ŷգ�
    !��̬������ˮ�����飨�£�
    REAL,ALLOCATABLE::h_new(:)

    !/////////////�Ǻ㶨������
    !δ֪���洢
    REAL,ALLOCATABLE::ZZ(:,:)
    REAL,ALLOCATABLE::QQ(:,:)
    REAL,ALLOCATABLE::Pool_Volume(:,:)

    !��һ����һֱ���洢
    REAL,ALLOCATABLE::ZZOLD(:,:)
    REAL,ALLOCATABLE::QQOLD(:,:)
    REAL,ALLOCATABLE::ELE_Volume(:,:)
    !׷��ϵ��
    REAL,ALLOCATABLE::SS(:)
    REAL,ALLOCATABLE::TT(:)



    Real::GateBottom,GateWidth,BeginBottom
    real::pollu_Lenth
    real::arrive_time
    real::t       !ʱ�����,��λ��min
    real::D  !��ȾԴ���ˮ��500�״�����
    real::S  !����õ���ȥȾ��ǰ�˾���ȾԴ����
    real::differ    !���
    real::D1,D2,w    !��������ó���D
    REAL::v,m,B,h,JJ,g   !ϵ��
    real,allocatable::v_eachTime(:)!ÿ�μ�����ٶ�
    real,allocatable::h_eachTime(:)!ÿ�μ����h

    allocate(v_eachTime(Boundary_cal.NUM))
    allocate(h_eachTime(Boundary_cal.NUM))
    !allocate(Opendegree(1,Boundary_cal.NUM))

    v=1.0
    m=1.0    !���ɱ���ȡ1
    B=26.5
    !h=6.25
    JJ=0.00001
    g=pollut_weight !��Ⱦ��������������Ϊ��10t����ˮ�����㵥λ:g,�������ֶ�������Ⱦ������
    D=pollut_ConveyDistance    !��λ:m

    !////////////////////////////////////////////////////
    NUMMember=Config_main.Assimilate_set.NUMMember
    m_coeff=Config_main.Assimilate_set.m_coeff
    PPQ=Config_main.Assimilate_set.PPQ
    VVQ=Config_main.Assimilate_set.VVQ
    Canal_NUM=Com_element_main.Canal_NUM
    c_element=Com_element_cal.Element_NUM

    !//////////////////////////////////
    Dynamic_Result_cal.DeltaT=Boundary_cal.Step
    Dynamic_Result_cal.NUM=Boundary_cal.NUM
    CALL Dynamic_Result_Auto_Update(Dynamic_Result_cal,Com_element_main,Boundary_cal,"NUM")
    ALLOCATE(h_new(c_element+1))
    !δ֪���洢
    ALLOCATE(ZZ(c_element+1,NUMMember))
    ALLOCATE(QQ(c_element+1,NUMMember))
    ALLOCATE(ZZOLD(c_element+1,NUMMember))
    ALLOCATE(QQOLD(c_element+1,NUMMember))
    ALLOCATE(SG_POSTION(Com_element_cal.Controlgate_NUM))
    ALLOCATE(DM_POSTION(Com_element_cal.Controlgate_NUM+2))

    ALLOCATE(Dynamic_Result_cal.Pool_Volume(Com_element_cal.Controlgate_NUM+1,Boundary_cal.NUM))
    ALLOCATE(ElE_Volume(c_element,NUMMember))
    ALLOCATE(Dynamic_Result_cal.Total_Volume(Boundary_cal.NUM))
    !///////ȷ��բ�����ڵ�λ��//////
    pos_sg=1
    DO j=1,c_element
        IF(Com_element_cal.Ele_Rel(j,2)==6)THEN
            SG_POSTION(pos_sg)=j
            pos_sg=pos_sg+1
        ENDIF
    END DO

    !///////ȷ��բ��+��β�������ڵ�λ��//////
    DM_sg=2
    DM_POSTION(1)=1
    DO j=2,c_element
        IF(Com_element_cal.Ele_Rel(j,2)==6)THEN
            DM_POSTION(DM_sg)=j
            DM_sg=DM_sg+1
        END IF
    END DO
    DM_POSTION(Com_element_cal.Controlgate_NUM+2)=c_element

    DO j=c_element+1,1,-1
        h_new(j)=Condition_0.h(c_element+2-j)
        QQ(j,1)=Condition_0.Q(c_element+2-j)
    END DO

    DO mm=1,NUMMember
        DO j=1,c_element,1
            IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
                ZZ(j,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j,3)).InInvertElev+h_new(j)
                ZZ(j+1,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j,3)).OutInvertElev+h_new(j+1)
            ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////������Ԫ��
                ZZ(j,mm)=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j,3)).InInvertElev+h_new(j)
                ZZ(j+1,mm)=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j,3)).OutInvertElev+h_new(j+1)
            ENDIF
        END DO
    END DO
    IF(Boundary_cal.Upsign==1)THEN
        ZZ(1,1)=Boundary_cal.UpZ(1)
    ENDIF

    !////////�ܵļ������Ϊ��Boundary_cal.NUM,�������ļ��м������
    si_n=1
    DO k=1,Boundary_cal.NUM
        DO mm=1,NUMMember
            DO i=1,c_element,1
                QQOLD(i,mm)=QQ(i,mm)
                if((Com_element_cal.Ele_Rel(i,2)==5).and.( ZZ(i,mm)<Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(i,3)).InInvertElev+0.1)) then
                    ZZOLD(i,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(i,3)).InInvertElev+0.1
                    write(*,*)"==============warning��׷�Ϸ���������==============="
                    write(*,*) "��",k,"�μ����ʱ��¶����,��ˮ�ֵΪ0.1��������"
                    write(*,*) "����������������¶�����⣬�ᵼ�½������������ԭ����δ��֪"
                    write(*,*) "���Գ��Խ����㲽�����̣����ߵ���һ�²���"
                    write(*,*)"================================================"
                else
                    ZZOLD(i,mm)=ZZ(i,mm)
                end if
            END DO
            !���һ�����浥������
            QQOLD(c_element+1,mm)=QQ(c_element+1,mm)
            if((Com_element_cal.Ele_Rel(c_element,2)==5).and.( ZZ(c_element,mm)<Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(c_element,3)).OutInvertElev+0.1)) then
                ZZOLD(c_element+1,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(c_element,3)).OutInvertElev+0.1
            else
                ZZOLD(c_element+1,mm)=ZZ(c_element+1,mm)
            end if
            !#TU-------------------------�޸ĵ��˽���
            !׷��ϵ��
            ALLOCATE(SS(c_element+1))
            ALLOCATE(TT(c_element+1))
            ALLOCATE(Condition_cal.PP(c_element+1))
            ALLOCATE(Condition_cal.VV(c_element+1))


            !***************����ϵ��************
            IF(Boundary_cal.Upsign==0)THEN        !���α߽������ֱ�ѡ��ˮλ������.0��ʾѡ�������߽磬1��ʾ����Ϊˮλ�߽�
                !/////////���α߽�ϵ��
                Condition_cal.PP(1)=Boundary_cal.UpQ(k)
                QQ(1,mm)=Boundary_cal.UpQ(k)
                Condition_cal.VV(1)=0.0
                !/////////�м�Ԫ��ϵ��
                DO j=1,c_element
                    !/////���ж������ٵ��ö�Ӧ�ĺ���ģ�飬������Ӧ��ϵ����
                    IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
                        Call CanalDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==4)THEN     !////��ˮ��Ԫ��
                        Call SideOutFlowDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_cal.OutDischarge(Com_element_cal.Ele_Rel(j,3),k),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==6)THEN     !////����բԪ��
                        Call ControlGateDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_temp.Opendegree(1,k),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                        firstgate_Index = j
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==7)THEN  !/////�����Ԫ��
                        Call AreachangeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////������Ԫ��
                        Call InvertedsiphonDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==11)THEN !/////�Ŷ�Ԫ��
                        Call BridgeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==12)THEN    !/////��վԪ��
                        Call PumpingstationDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),C_Bladeangle_cl(Com_element_cal.Ele_Rel(j,3)).Each_PS(:,k),Com_element_cal.PumpST_cl(Com_element_cal.Ele_Rel(j,3)).UNIT_NUM,QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ENDIF
                END DO
                !///////////���α߽�����
                IF(Boundary_cal.Downsign==1)THEN   !///////////ˮλ�߽�1
                    ZZ(c_element+1,mm)=Boundary_cal.DownZ(k)
                ELSEIF(Boundary_cal.Downsign==2)THEN   !///////////ˮλ~������ϵ2    �˴������α߽��ʹ�õ�DownZ�����ܻ�������
                    PPQ=m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))-1.0/(2.0*m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))
                    VVQ=-1.0/(2.0*m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev)))
                    ZZ(c_element+1,mm)=(Condition_cal.PP(c_element+1)-PPQ)/(Condition_cal.VV(c_element+1)-VVQ)
                ENDIF

                !*******�ش����ˮ��������********��
                DO j=c_element,1,-1
                    ZZ(j,mm)=SS(j+1)-TT(j+1)*ZZ(j+1,mm)
                    QQ(j+1,mm)=Condition_cal.PP(j+1)-Condition_cal.VV(j+1)*ZZ(j+1,mm)
                END DO
            ELSE    !����Ϊˮλ�߽�Boundary_cal.Upsign==1
                !/////////���α߽�ϵ��
                Condition_cal.PP(1)=Boundary_cal.UpZ(k)
                ZZ(1,mm)=Boundary_cal.UpZ(k)
                Condition_cal.VV(1)=0.0
                !/////////�м�Ԫ��ϵ��
                DO j=1,c_element
                    !/////���ж������ٵ��ö�Ӧ�ĺ���ģ�飬������Ӧ��ϵ����
                    IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
                        Call CanalDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==4)THEN     !////��ˮ��Ԫ��
                        Call SideOutFlowDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_cal.OutDischarge(Com_element_cal.Ele_Rel(j,3),k),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==6)THEN     !////����բԪ��
                        Call ControlGateDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_temp.Opendegree(1,k),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                        firstgate_Index = j
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==7)THEN  !/////�����Ԫ��
                        Call AreachangeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////������Ԫ��
                        Call InvertedsiphonDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==11)THEN !/////�Ŷ�Ԫ��
                        Call BridgeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==12)THEN    !/////��վԪ��
                        Call PumpingstationDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),C_Bladeangle_cl(Com_element_cal.Ele_Rel(j,3)).Each_PS(:,k),Com_element_cal.PumpST_cl(Com_element_cal.Ele_Rel(j,3)).UNIT_NUM,QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ENDIF
                END DO
                !///////////���α߽�����
                IF(Boundary_cal.Downsign==0)THEN   !///////////�����߽�0
                    QQ(c_element+1,mm)=Boundary_cal.DownQ(k)
                ELSEIF(Boundary_cal.Downsign==1)THEN   !///////////ˮλ�߽�1
                    QQ(c_element+1,mm)=(Condition_cal.PP(c_element+1)-Boundary_cal.DownZ(k))/Condition_cal.VV(c_element+1)
                ELSEIF(Boundary_cal.Downsign==2)THEN   !///////////ˮλ~������ϵ�߽�2
                    PPQ=m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))-1.0/(2.0*m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))
                    VVQ=-1.0/(2.0*m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev)))
                    QQ(c_element+1,mm)=(PPQ-VVQ*Condition_cal.PP(c_element+1))/(1-VVQ*Condition_cal.VV(c_element+1))
                ENDIF
                !*******�ش����ˮ��������********��
                DO j=c_element,1,-1
                    QQ(j,mm)=SS(j+1)-TT(j+1)*QQ(j+1,mm)
                    ZZ(j+1,mm)=Condition_cal.PP(j+1)-Condition_cal.VV(j+1)*QQ(j+1,mm)
                END DO

            END IF

            IF(Output_sign==0)THEN
                DO i=1,Com_element_cal.Controlgate_NUM
                    IF(MOD((k-1)*Boundary_cal.Step,Out_DeltaT)==0)THEN
                        Dynamic_Result_cal.Zup(i,k)=ZZ(sg_postion(i),mm)
                        Dynamic_Result_cal.Zdown(i,k)=ZZ(sg_postion(i)+1,mm)
                        Dynamic_Result_cal.Q(i,k)=QQ(sg_postion(i),mm)
                    END IF
                END DO
            ELSEIF(Output_sign==1)THEN
                DO i=1,Com_element_main.Element_NUM-Com_element_main.Bridge_NUM
                    IF(MOD((k-1)*Boundary_cal.Step,Out_DeltaT)==0)THEN
                        Dynamic_Result_cal.ZZ(i,k)=ZZ(Com_element_cal.ELE_POSTION(i),mm)
                        Dynamic_Result_cal.QQ(i,k)=QQ(Com_element_cal.ELE_POSTION(i),mm)
                    END IF
                END DO
            END IF

            Dynamic_Result_cal.Zbd(1,k)=ZZ(c_element+1,mm)
            Dynamic_Result_cal.Qbd(1,k)=QQ(c_element+1,mm)

            !***���������ܺ�*******
            DO i=1,Com_element_cal.Controlgate_NUM+1
                Dynamic_Result_cal.Pool_Volume(i,k)=0.0
                DO j=DM_POSTION(i),DM_POSTION(i+1)-1
                    Dynamic_Result_cal.Pool_Volume(i,k)=Dynamic_Result_cal.Pool_Volume(i,k)+ElE_Volume(j,mm)
                END DO
            END DO
            !Dynamic_Result_cal.Pool_Volume(Controlgate_NUM+1,k)=Dynamic_Result_cal.Pool_Volume(Controlgate_NUM+1,k)+ElE_Volume(c_element,mm)

            !***�����ܺ�*******
            Dynamic_Result_cal.Total_Volume(k)=0.0
            DO i=1,Com_element_cal.Controlgate_NUM+1
                Dynamic_Result_cal.Total_Volume(k)=Dynamic_Result_cal.Total_Volume(k)+Dynamic_Result_cal.Pool_Volume(i,k)
            END DO

            GateBottom=Com_element_cal.Cogate_cl(1).InvertElev
            GateWidth=Com_element_cal.Cogate_cl(1).ParallelNum*Com_element_cal.Cogate_cl(1).SingleWidth
            BeginBottom=Com_element_cal.Ca_cl(1).InInvertElev
            !��һ���ĺӽ���բǰ����
            v_eachTime(k)=1.2*QQ(firstgate_Index,mm)/(ZZ(firstgate_Index,mm)-GateBottom)/GateWidth!ϵ������ƽ����ˮ�������½�
            !��һ�µ��բ��ˮλ���ĺ�բǰˮλ��ƽ��ֵ
            h_eachTime(k)=((ZZ(1,mm)-BeginBottom)+(ZZ(firstgate_Index,mm)-GateBottom))/2
            !!47Ϊ�ĺ�բǰ�������Ϊ76
            !v_eachTime(k)=1.2*QQ(76,mm)/(ZZ(76,mm)-138.53)/26!����բǰ���٣�27�׿�139.25�׸̣߳����������ˣ���΢�˸�ϵ��ƽ����ˮ�����ٽ���
            !!��һ�µ��բ��ˮλ�������բǰˮλ��ƽ��ֵ
            !h_eachTime(k)=((ZZ(1,mm)-140.7)+(ZZ(76,mm)-138.53))/2


            DEALLOCATE(SS)
            DEALLOCATE(TT)
            DEALLOCATE(Condition_cal.PP)
            DEALLOCATE(Condition_cal.VV)
        END DO
    END DO

    !////////////////////////////////////////////////////////////////////////////////////////////////

    v=sum(v_eachTime(:))/Boundary_cal.NUM
    h=sum(h_eachTime(:))/Boundary_cal.NUM
    !��Ⱦ������
    DO t=1,1000,0.01
        D1=60*v*t
        D2=(m*0.011*(v**2*B**2))/(h*SQRT(9.8*h*JJ))
        w=(12+log(g/10))*SQRT(2*D2)*t**0.455
        S=D1+W/2
        differ=D-S
        if(differ>10)then
            cycle
        else if(differ<10)then
            exit
        end if
    end do
    pollu_Lenth=w
    arrive_time=t
    DEALLOCATE(v_eachTime)
    DEALLOCATE(h_eachTime)


    !///////�ͷſռ�
    DEALLOCATE(h_new)
    !DEALLOCATE(ZZ)
    !DEALLOCATE(QQ)
    DEALLOCATE(ZZOLD)
    DEALLOCATE(QQOLD)
    DEALLOCATE(SG_POSTION)
    END















    !���Ĳ���
    SUBROUTINE Dynamic_Cal2(Boundary_temp,pollu_Lenth,arrive_time,Endtime)  !�Ǻ㶨������
    USE GLOBAL_VARIABLE

    IMPLICIT NONE

    TYPE(Boundary)::Boundary_temp


    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_main,Com_element_cal
    !TYPE(Config)::Config_main

    TYPE(Dynamic_Result)::Dynamic_Result_cal
    !TYPE(Condition)::Condition_0,Condition_cal
    TYPE(Condition)::Condition_cal
    INTEGER::NUMMember
    INTEGER,ALLOCATABLE::DM_POSTION(:)  !����բ+��β����λ��
    INTEGER,ALLOCATABLE::SG_POSTION(:)  !����բλ��
    REAL::m_coeff    !ĩ��բ�Ź�բ����ϵ������
    REAL::PPQ, VVQ  !ˮλ~������ϵ����ϵ��
    INTEGER::Canal_NUM
    INTEGER::c_element    !��������ռ䲽����������Ԫ���ĸ���

    INTEGER::i,j,k,mm,kk,LLL
    real::c11,c22,c21,hh
    REAL::zzaverage,QQaverage
    INTEGER::si_n   !�Ż��洢ʱ��ѭ������
    REAL::random    !�����
    DOUBLE PRECISION r0
    Integer::pos_sg,DM_sg !��¼����բλ��
    INTEGER::ele_sg    !��¼��Ԫ����λ�ã������Ŷգ�
    !��̬������ˮ�����飨�£�
    REAL,ALLOCATABLE::h_new(:)

    !/////////////�Ǻ㶨������
    !δ֪���洢
    REAL,ALLOCATABLE::ZZ(:,:)
    REAL,ALLOCATABLE::QQ(:,:)
    REAL,ALLOCATABLE::Pool_Volume(:,:)

    !��һ����һֱ���洢
    REAL,ALLOCATABLE::ZZOLD(:,:)
    REAL,ALLOCATABLE::QQOLD(:,:)
    REAL,ALLOCATABLE::ELE_Volume(:,:)
    !׷��ϵ��
    REAL,ALLOCATABLE::SS(:)
    REAL,ALLOCATABLE::TT(:)



    Real::GateBottom,GateWidth,BeginBottom

    real::pollu_Lenth
    real::arrive_time
    real::t       !ʱ�����,��λ��min
    real::D  !��ȾԴ���ˮ��500�״�����
    real::S  !����õ���ȥȾ��ǰ�˾���ȾԴ����
    real::differ    !���
    real::D1,D2,w    !��������ó���D
    REAL::v,m,B,h,JJ,g   !ϵ��
    real,allocatable::v_eachTime(:)!ÿ�μ�����ٶ�
    real,allocatable::h_eachTime(:)!ÿ�μ����h

    real,allocatable :: Q_char(:),Z_char(:)
    CHARACTER(LEN=100):: temp
    real::temp2
    Integer::Endtime,isovergate,mindex
    Real,allocatable::gatedischarge(:)
    Integer::canal_f
    Integer::invert_f


    allocate(v_eachTime(Boundary_2cal.NUM))
    allocate(h_eachTime(Boundary_2cal.NUM))
    !allocate(Opendegree(1,Boundary_2cal.NUM))
    v=1.0
    m=1.0    !���ɱ���ȡ1
    B=26.5
    h=6.25
    JJ=0.00001
    g=pollut_weight !��Ⱦ��������������Ϊ��10t����ˮ�����㵥λ:g,�������ֶ�������Ⱦ������
    D=pollut_ConveyDistance    !��λ:m



    !////////////////////////////////////////////////////
    NUMMember=Config_main.Assimilate_set.NUMMember
    m_coeff=Config_main.Assimilate_set.m_coeff
    PPQ=Config_main.Assimilate_set.PPQ
    VVQ=Config_main.Assimilate_set.VVQ
    Canal_NUM=Com_element_main.Canal_NUM
    c_element=Com_element_cal.Element_NUM

    !//////////////////////////////////
    Dynamic_Result_cal.DeltaT=Boundary_2cal.Step
    Dynamic_Result_cal.NUM=Boundary_2cal.NUM
    CALL Dynamic_Result_Auto_Update(Dynamic_Result_cal,Com_element_main,Boundary_2cal,"NUM")
    ALLOCATE(h_new(c_element+1))
    !δ֪���洢
    ALLOCATE(ZZ(c_element+1,NUMMember))
    ALLOCATE(QQ(c_element+1,NUMMember))
    ALLOCATE(ZZOLD(c_element+1,NUMMember))
    ALLOCATE(QQOLD(c_element+1,NUMMember))
    ALLOCATE(SG_POSTION(Com_element_cal.Controlgate_NUM))
    ALLOCATE(DM_POSTION(Com_element_cal.Controlgate_NUM+2))

    ALLOCATE(Dynamic_Result_cal.Pool_Volume(Com_element_cal.Controlgate_NUM+1,Boundary_2cal.NUM))
    ALLOCATE(ElE_Volume(c_element,NUMMember))
    ALLOCATE(Dynamic_Result_cal.Total_Volume(Boundary_2cal.NUM))
    !///////ȷ��բ�����ڵ�λ��//////
    pos_sg=1
    DO j=1,c_element
        IF(Com_element_cal.Ele_Rel(j,2)==6)THEN
            SG_POSTION(pos_sg)=j
            pos_sg=pos_sg+1
        ENDIF
    END DO

    !///////ȷ��բ��+��β�������ڵ�λ��//////
    DM_sg=2
    DM_POSTION(1)=1
    DO j=2,c_element
        IF(Com_element_cal.Ele_Rel(j,2)==6)THEN
            DM_POSTION(DM_sg)=j
            DM_sg=DM_sg+1
        END IF
    END DO
    DM_POSTION(Com_element_cal.Controlgate_NUM+2)=c_element


    DO j=c_element+1,1,-1
        h_new(j)=Condition_0.h(c_element+2-j)
        QQ(j,1)=Condition_0.Q(c_element+2-j)
    END DO

    DO mm=1,NUMMember
        DO j=1,c_element,1
            IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
                ZZ(j,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j,3)).InInvertElev+h_new(j)
                ZZ(j+1,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j,3)).OutInvertElev+h_new(j+1)
            ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////������Ԫ��
                ZZ(j,mm)=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j,3)).InInvertElev+h_new(j)
                ZZ(j+1,mm)=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j,3)).OutInvertElev+h_new(j+1)
            ENDIF
        END DO
    END DO
    IF(Boundary_2cal.Upsign==1)THEN
        ZZ(1,1)=Boundary_2cal.UpZ(1)
    ENDIF

    !////////�ܵļ������Ϊ��Boundary_2cal.NUM,�������ļ��м������
    si_n=1
    DO k=1,Boundary_2cal.NUM
        DO mm=1,NUMMember

            DO i=1,c_element,1
                QQOLD(i,mm)=QQ(i,mm)
                if((Com_element_cal.Ele_Rel(i,2)==5).and.( ZZ(i,mm)<Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(i,3)).InInvertElev+0.1)) then
                    ZZOLD(i,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(i,3)).InInvertElev+0.1
                    write(*,*)"==============warning��׷�Ϸ���������==============="
                    write(*,*) "��",k,"�μ����ʱ��¶����,��ˮ�ֵΪ0.1��������"
                    write(*,*) "����������������¶�����⣬�ᵼ�½������������ԭ����δ��֪"
                    write(*,*) "���Գ��Խ����㲽�����̣����ߵ���һ�²���"
                    write(*,*)"================================================"
                else
                    ZZOLD(i,mm)=ZZ(i,mm)
                end if
            END DO


            !���һ�����浥������
            QQOLD(c_element+1,mm)=QQ(c_element+1,mm)
            if((Com_element_cal.Ele_Rel(c_element,2)==5).and.( ZZ(c_element,mm)<Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(c_element,3)).OutInvertElev+0.1)) then
                ZZOLD(c_element+1,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(c_element,3)).OutInvertElev+0.1
            else
                ZZOLD(c_element+1,mm)=ZZ(c_element+1,mm)
            end if







            !#TU-------------------------�޸ĵ��˽���
            !׷��ϵ��
            ALLOCATE(SS(c_element+1))
            ALLOCATE(TT(c_element+1))
            ALLOCATE(Condition_cal.PP(c_element+1))
            ALLOCATE(Condition_cal.VV(c_element+1))







            isovergate = 0
            !***************����ϵ��************
            IF(Boundary_2cal.Upsign==0)THEN        !���α߽������ֱ�ѡ��ˮλ������.0��ʾѡ�������߽磬1��ʾ����Ϊˮλ�߽�
                !/////////���α߽�ϵ��
                Condition_cal.PP(1)=Boundary_2cal.UpQ(k)
                QQ(1,mm)=Boundary_2cal.UpQ(k)
                Condition_cal.VV(1)=0.0
                !/////////�м�Ԫ��ϵ��
                DO j=1,c_element
                    if((isovergate == 1).and.(k > Endtime))Then
                        QQOLD(j:,m)=150
                        ZZOLD(j,m)=ZZOLD(j,m)+0.6
                        isovergate = 0
                    ENDIF

                    !/////���ж������ٵ��ö�Ӧ�ĺ���ģ�飬������Ӧ��ϵ����
                    IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
                        Call CanalDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==4)THEN     !////��ˮ��Ԫ��
                        Call SideOutFlowDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_2cal.OutDischarge(Com_element_cal.Ele_Rel(j,3),k),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==6)THEN     !////����բԪ��
                        Call ControlGateDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_temp.Opendegree(1,k),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                        isovergate = 1
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==7)THEN  !/////�����Ԫ��
                        Call AreachangeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////������Ԫ��
                        Call InvertedsiphonDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==11)THEN !/////�Ŷ�Ԫ��
                        Call BridgeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==12)THEN    !/////��վԪ��
                        Call PumpingstationDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),C_Bladeangle_cl(Com_element_cal.Ele_Rel(j,3)).Each_PS(:,k),Com_element_cal.PumpST_cl(Com_element_cal.Ele_Rel(j,3)).UNIT_NUM,QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ENDIF
                END DO
                !///////////���α߽�����
                IF(Boundary_2cal.Downsign==1)THEN   !///////////ˮλ�߽�1
                    ZZ(c_element+1,mm)=Boundary_2cal.DownZ(k)
                ELSEIF(Boundary_2cal.Downsign==2)THEN   !///////////ˮλ~������ϵ2    �˴������α߽��ʹ�õ�DownZ�����ܻ�������
                    PPQ=m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))-1.0/(2.0*m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))
                    VVQ=-1.0/(2.0*m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev)))
                    ZZ(c_element+1,mm)=(Condition_cal.PP(c_element+1)-PPQ)/(Condition_cal.VV(c_element+1)-VVQ)
                ENDIF

                !*******�ش����ˮ��������********��
                DO j=c_element,1,-1
                    ZZ(j,mm)=SS(j+1)-TT(j+1)*ZZ(j+1,mm)
                    QQ(j+1,mm)=Condition_cal.PP(j+1)-Condition_cal.VV(j+1)*ZZ(j+1,mm)
                END DO
            ELSE    !����Ϊˮλ�߽�Boundary_2cal.Upsign==1
                !/////////���α߽�ϵ��
                Condition_cal.PP(1)=Boundary_2cal.UpZ(k)
                ZZ(1,mm)=Boundary_2cal.UpZ(k)
                Condition_cal.VV(1)=0.0
                !/////////�м�Ԫ��ϵ��
                DO j=1,c_element
                    !/////���ж������ٵ��ö�Ӧ�ĺ���ģ�飬������Ӧ��ϵ����
                    IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////����Ԫ��
                        Call CanalDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==4)THEN     !////��ˮ��Ԫ��
                        Call SideOutFlowDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_2cal.OutDischarge(Com_element_cal.Ele_Rel(j,3),k),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==6)THEN     !////����բԪ��
                        Call ControlGateDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_temp.Opendegree(1,k),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                        isovergate = 1
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==7)THEN  !/////�����Ԫ��
                        Call AreachangeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////������Ԫ��
                        Call InvertedsiphonDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==11)THEN !/////�Ŷ�Ԫ��
                        Call BridgeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==12)THEN    !/////��վԪ��
                        Call PumpingstationDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),C_Bladeangle_cl(Com_element_cal.Ele_Rel(j,3)).Each_PS(:,k),Com_element_cal.PumpST_cl(Com_element_cal.Ele_Rel(j,3)).UNIT_NUM,QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ENDIF
                END DO
                !///////////���α߽�����
                IF(Boundary_2cal.Downsign==0)THEN   !///////////�����߽�0
                    QQ(c_element+1,mm)=Boundary_2cal.DownQ(k)
                ELSEIF(Boundary_2cal.Downsign==1)THEN   !///////////ˮλ�߽�1
                    QQ(c_element+1,mm)=(Condition_cal.PP(c_element+1)-Boundary_2cal.DownZ(k))/Condition_cal.VV(c_element+1)
                ELSEIF(Boundary_2cal.Downsign==2)THEN   !///////////ˮλ~������ϵ�߽�2
                    PPQ=m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))-1.0/(2.0*m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))
                    VVQ=-1.0/(2.0*m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev)))
                    QQ(c_element+1,mm)=(PPQ-VVQ*Condition_cal.PP(c_element+1))/(1-VVQ*Condition_cal.VV(c_element+1))
                ENDIF
                !*******�ش����ˮ��������********��
                DO j=c_element,1,-1
                    QQ(j,mm)=SS(j+1)-TT(j+1)*QQ(j+1,mm)
                    ZZ(j+1,mm)=Condition_cal.PP(j+1)-Condition_cal.VV(j+1)*QQ(j+1,mm)
                END DO

            END IF

            IF(Output_sign==0)THEN
                DO i=1,Com_element_cal.Controlgate_NUM
                    IF(MOD((k-1)*Boundary_2cal.Step,Out_DeltaT)==0)THEN
                        Dynamic_Result_cal.Zup(i,k)=ZZ(sg_postion(i),mm)
                        Dynamic_Result_cal.Zdown(i,k)=ZZ(sg_postion(i)+1,mm)
                        Dynamic_Result_cal.Q(i,k)=QQ(sg_postion(i),mm)
                    END IF
                END DO
            ELSEIF(Output_sign==1)THEN
                DO i=1,Com_element_main.Element_NUM-Com_element_main.Bridge_NUM
                    IF(MOD((k-1)*Boundary_2cal.Step,Out_DeltaT)==0)THEN
                        Dynamic_Result_cal.ZZ(i,k)=ZZ(Com_element_cal.ELE_POSTION(i),mm)
                        Dynamic_Result_cal.QQ(i,k)=QQ(Com_element_cal.ELE_POSTION(i),mm)
                    END IF
                END DO
            END IF

            Dynamic_Result_cal.Zbd(1,k)=ZZ(c_element+1,mm)
            Dynamic_Result_cal.Qbd(1,k)=QQ(c_element+1,mm)

            !***���������ܺ�*******
            DO i=1,Com_element_cal.Controlgate_NUM+1
                Dynamic_Result_cal.Pool_Volume(i,k)=0.0
                DO j=DM_POSTION(i),DM_POSTION(i+1)-1
                    Dynamic_Result_cal.Pool_Volume(i,k)=Dynamic_Result_cal.Pool_Volume(i,k)+ElE_Volume(j,mm)
                END DO
            END DO
            !Dynamic_Result_cal.Pool_Volume(Controlgate_NUM+1,k)=Dynamic_Result_cal.Pool_Volume(Controlgate_NUM+1,k)+ElE_Volume(c_element,mm)

            !***�����ܺ�*******
            Dynamic_Result_cal.Total_Volume(k)=0.0
            DO i=1,Com_element_cal.Controlgate_NUM+1
                Dynamic_Result_cal.Total_Volume(k)=Dynamic_Result_cal.Total_Volume(k)+Dynamic_Result_cal.Pool_Volume(i,k)
            END DO



            GateBottom=Com_element_cal.Cogate_cl(1).InvertElev
            GateWidth=Com_element_cal.Cogate_cl(1).ParallelNum*Com_element_cal.Cogate_cl(1).SingleWidth
            BeginBottom=Com_element_cal.Ca_cl(1).InInvertElev
            v_eachTime(k)=1.2*QQ(firstgate_Index,mm)/(ZZ(firstgate_Index,mm)-GateBottom)/GateWidth!ϵ������ƽ����ˮ�������½�
            h_eachTime(k)=((ZZ(1,mm)-BeginBottom)+(ZZ(firstgate_Index,mm)-GateBottom))/2




            !*******��մ洢�ռ�WRITE(24,*)
            DEALLOCATE(SS)
            DEALLOCATE(TT)
            DEALLOCATE(Condition_cal.PP)
            DEALLOCATE(Condition_cal.VV)
        END DO
    END DO

    !////////////////////////////////////////////////////////////////////////////////////////////////


    v=sum(v_eachTime(:))/Boundary_2cal.NUM
    h=sum(h_eachTime(:))/Boundary_2cal.NUM
    !��Ⱦ������
    DO t=1,1000,0.01
        D1=60*v*t
        D2=(m*0.011*(v**2*B**2))/(h*SQRT(9.8*h*JJ))
        w=(12+log(g/10))*SQRT(2*D2)*t**0.455
        S=D1+W/2
        differ=D-S
        if(differ>10)then
            cycle
        else if(differ<10)then
            exit
        end if
    end do
    pollu_Lenth=w
    arrive_time=t
    DEALLOCATE(v_eachTime)
    DEALLOCATE(h_eachTime)
    !///////�ͷſռ�
    DEALLOCATE(h_new)
    !DEALLOCATE(ZZ)
    !DEALLOCATE(QQ)
    DEALLOCATE(ZZOLD)
    DEALLOCATE(QQOLD)
    DEALLOCATE(SG_POSTION)


    !******************����Ǻ㶨������ˮ��������******************
    ALLOCATE (Q_char(firstGateIndex))
    ALLOCATE (Z_char(firstGateIndex))
    canal_f=0
    invert_f=0
    DO i = 1,firstGateIndex
        if(Com_element_main.ele_rel(I,2)==5)then
            canal_f=canal_f+1
            Q_char(i)= Com_element_main.Ca_cl(canal_f).INcoord
            Q_char(i+1)= Com_element_main.Ca_cl(canal_f).Outcoord
        ENDIF
        if(Com_element_main.ele_rel(I,2)==10)then
            invert_f=invert_f+1
            Q_char(i)= Com_element_main.insi_cl(invert_f).INcoord
            Q_char(i+1)= Com_element_main.insi_cl(invert_f).Outcoord
        ENDIF
    END DO

    
    !!��������
    OPEN(UNIT=28,FILE="OUTPUTFILE/DynamicRESULT_Q.txt")
    WRITE(28,'(*(G0.5,:,"   "))')'ʱ��(����)\׮�ţ�km��',(Q_char(i),i=1,firstGateIndex)




    temp2=Dynamic_Result_cal.ZZ(10,1)
    DO k=1,Dynamic_Result_cal.NUM,1
        DO i=1,firstGateIndex
            IF ((Dynamic_Result_cal.QQ(i,k)<0.5).and.(Dynamic_Result_cal.QQ(i,k)>-0.5))Then
                Dynamic_Result_cal.QQ(i,k)=0
            ENDIF
        ENDDO

    END DO


    allocate(gatedischarge(size(Boundary_2in.Opendegree,2)))

    mindex = 0
    DO k=1,Dynamic_Result_cal.NUM,1
        IF(MOD((k-1)*Dynamic_Result_cal.DeltaT,Out_DeltaT)==0)THEN

            !WRITE(28,'(*(G0.5,:,","))')(k-1)*Dynamic_Result_cal.DeltaT/REAL(Out_DeltaT),(Dynamic_Result_cal.ZZ(i,k),Dynamic_Result_cal.QQ(i,k),i=1,Com_element_main.Element_NUM-Com_element_main.Bridge_NUM),Dynamic_Result_cal.Zbd(1,k),Dynamic_Result_cal.Qbd(1,k)!����
            WRITE(28,'(*(G0.5,:,"   "))')(k-1)*Dynamic_Result_cal.DeltaT/REAL(Out_DeltaT)*Boundary_in.Step/60.0,(Dynamic_Result_cal.QQ(i,k),i=1,firstGateIndex)
            mindex=mindex+1
            gatedischarge(mindex)=Dynamic_Result_cal.QQ(firstGateIndex,k)

            if (temp2>Dynamic_Result_cal.ZZ(10,k)) temp2=Dynamic_Result_cal.ZZ(10,k)
        END IF
    END DO
    CLOSE(28)


    !!ˮλ����
    OPEN(UNIT=28,FILE="OUTPUTFILE/DynamicRESULT_Z.txt")
    if(Output_sign==1)then
        !WRITE(28,'(*(G0.5,:,"  "))')'ʱ������',(Z_char(i),Q_char(i),i=1,Com_element_main.Element_NUM-Com_element_main.Bridge_NUM+1)
        WRITE(28,'(*(G0.5,:,"   "))')'ʱ��(����)\׮�ţ�km��',(Q_char(i),i=1,firstGateIndex)
    end if
    DO k=1,Dynamic_Result_cal.NUM,1
        IF(MOD((k-1)*Dynamic_Result_cal.DeltaT,Out_DeltaT)==0)THEN
            WRITE(28,'(*(G0.5,:,"   "))')(k-1)*Dynamic_Result_cal.DeltaT/REAL(Out_DeltaT)*Boundary_in.Step/60.0,(Dynamic_Result_cal.ZZ(i,k),i=1,firstGateIndex)
        END IF
    END DO
    CLOSE(28)


    !����Ż���Ľ���բ����
    open(unit=71,file="outputfile/selected-operating.txt")
    write(71,*)"ʱ��(min) �½���բ����(mm)  �½���բQ  ��ˮբQ(m3/s) �ܳ���Q   �Ͻ���բQ  "
    do i=1,size(Boundary_2in.Opendegree,2)
        WRITE(71,'(*(G0.5,:,"   "))') (i-1)*Boundary_2in.Step/60,int(Boundary_2in.Opendegree(1,i)*1000),gatedischarge(i),Boundary_2in.OutDischarge(outflow_index,i),gatedischarge(i)+Boundary_2in.OutDischarge(outflow_index,i),Boundary_2in.UPQ(i)
    end do
    close(71)

    !OPEN(UNIT=29,FILE="OUTPUTFILE/out-ini.txt")
    !WRITE(29,'(*(G0.5,:,"   "))')(ZZ(i,1),i=1,53)
    !WRITE(29,'(*(G0.5,:,"   "))')(QQ(i,1),i=1,53)
    !CLOSE(29)

    !write(*,*) '��ˮբǰ���ˮλΪ:', temp2-137.642


    END
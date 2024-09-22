    SUBROUTINE Inputfile        !�����ļ�
    !**************************
    !******** ���¼�¼ ********
    !**************************
    !�Ż������ļ���� by Vincent Zhang 23/5/2022

    USE GLOBAL_VARIABLE

    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_in,Boundary_cal
    !TYPE(Com_element)::Com_element_main
    !TYPE(Config)::Config_main

    REAL::Up_tmp,Down_tmp

    INTEGER::i,j,k,PS_Co_NUM
    INTEGER::irr_num!#tuÿһ�������Ĳ�����������ͳ��
    INTEGER::section_num!#TU��Ҫ��ȡ�Ķ�����
    INTEGER::irrseries_num!#tu ������ ���������� ���������ܶ�β�����������β�����һ��
    INTEGER::irr_num_mat(200)!#tu
    INTEGER,allocatable::CDS(:)!#��ȡ����ʱ�ĸ̶߳���
    CHARACTER(len=12),allocatable::SectionName(:)!#tu��������
    REAL::GCD_temp(2,3000,200)!#һά��ά�Ǹ̵߳㣬��ά�Ƕ������
    INTEGER::kk,temp!#tu
    INTEGER::sum_of_irr!#tu��������������

    !**************************
    !******** ���б��� ********
    !**************************
    Config_main.Coeff_set.Theta=0.6
    Config_main.Coeff_set.Phi=0.5
    Config_main.Coeff_set.Psi=0.6
    Config_main.Coeff_set.PsiR=0.5
    Config_main.Coeff_set.Beta=1.0
    Config_main.Coeff_set.GRAV=9.81
    !RoughnessCoeff=0.02
    Config_main.Assimilate_set.NUMMember=1
    Config_main.Assimilate_set.m_coeff=0.35



    !/////�����ļ�
    OPEN(UNIT=18,FILE="INPUTFILE/Config.txt")
    READ(18,*) !����
    READ(18,*) !��Ⱦ������
    READ(18,*) pollut_weight
    READ(18,*) !�����Ͻ���բ����(m)
    READ(18,*) pollut_location
    READ(18,*) !��ˮբ���
    READ(18,*) outflow_index
    READ(18,*) !��ˮբ�������
    READ(18,*) outflow_maxQ
    READ(18,*) !��ˮբ�������
    READ(18,*) outflow_changelimit
    READ(18,*) !�Ͻ���բ������ƣ�������
    READ(18,*) inGate_changelimit
    READ(18,*) !����ʱ��Ȩ��
    READ(18,*) arrivetime_weight
    READ(18,*) !��Ⱥ����
    READ(18,*) populations
    READ(18,*) !��������
    READ(18,*) generations
    CLOSE(18)



    !**************************
    !******** �߽����� ********
    !**************************
    Boundary_in.Step=300
    Boundary_cal.Step=30
    Boundary_in.Start_Time=0
    Boundary_in.End_Time=24000!��ʱ�ֶ��̶�
    Out_DeltaT = Boundary_in.Step
    Boundary_2cal.Step=Boundary_cal.Step
    !Boundary_in.End_Time=Boundary_in.Start_Time+(NP/1)*Boundary_in.Step!�����2�����˽���բ������Ŀǰֻ���ֶ�������
    CALL Boundary_Auto_Update(Boundary_in,"NUM")
    !READ(11,*)
    !READ(11,*)Boundary_in.Upsign,Boundary_in.Downsign   !�߽��������Ʒ���0��ʾ�����߽磬1��ʾˮλ�߽磬2��ʾˮλ������ϵ�߽磨2���ܷ��������Σ�
    Boundary_in.Upsign=0
    Boundary_in.Downsign=1
    Output_sign=1
    !READ(11,*)
    !READ(11,*)Output_sign            !���������Ʒ���0-ֻ�������բ�������1-�����Ԫ��������
    !READ(11,*)
    !DO i=1,Boundary_in.NUM !���߽�����
    !    READ(11,*)Up_tmp,Down_tmp
    !    IF(Boundary_in.Upsign==0)THEN
    !        Boundary_in.UpQ(i)=Up_tmp
    !    ELSE IF(Boundary_in.Upsign==1)THEN
    !        Boundary_in.UpZ(i)=Up_tmp
    !    END IF
    !    IF(Boundary_in.Downsign==0.or.Boundary_in.Downsign==2)THEN
    !        Boundary_in.DownQ(i)=Down_tmp
    !    ELSE IF(Boundary_in.Downsign==1)THEN
    !        Boundary_in.DownZ(i)=Down_tmp
    !    END IF
    !END DO

    !//////����Ԫ��������Ϣ�����桢��̡����ȣ�
    OPEN(UNIT=12,FILE="INPUTFILE/COM-ELEMENT.txt")
    READ(12,*) !����Ԫ�������������໥���ӹ�ϵ�ļ�
    READ(12,*)!================================================================================
    READ(12,*) !Ԫ������
    READ(12,*)Com_element_main.Element_NUM

    READ(12,*) !�����͸���
    READ(12,*)
    READ(12,*)Com_element_main.Canal_NUM,Com_element_main.Areachange_NUM,Com_element_main.Invertedsiphon_NUM,&
        Com_element_main.Controlgate_NUM,Com_element_main.SideOutFlow_NUM,Com_element_main.Bridge_NUM,Com_element_main.PumpST_NUM
    CALL Com_element_Auto_Update(Com_element_main,"NUM")
    Boundary_in.Gate_NUM=Com_element_main.Controlgate_NUM
    CALL Boundary_Auto_Update(Boundary_in,"Gate_NUM")
    Boundary_in.PumpST_NUM=Com_element_main.PumpST_NUM
    CALL Boundary_Auto_Update(Boundary_in,"PumpST_NUM")
    Boundary_in.SideOutFlow_NUM=Com_element_main.SideOutFlow_NUM
    CALL Boundary_Auto_Update(Boundary_in,"SideOutFlow_NUM")
    READ(12,*)!================================================================================
    READ(12,*)

    DO i=1,Com_element_main.Element_NUM
        READ(12,*)(Com_element_main.Ele_Rel(i,j),j=1,3)
        if(Com_element_main.Ele_Rel(i,2)==6)Then
            firstGateIndex = Com_element_main.Ele_Rel(i,1)
        ENDIF
    END DO
    READ(12,*)!================================================================================
    !///������������-Ԫ������5
    READ(12,*)
    DO i=1,Com_element_main.Canal_NUM
        READ(12,*)Com_element_main.Ca_cl(i).Serial_NUM,&
            Com_element_main.Ca_cl(i).Name,&
            Com_element_main.Ca_cl(i).Incoord,&
            Com_element_main.Ca_cl(i).Outcoord,&
            Com_element_main.Ca_cl(i).BottomWidth, &              !wjb0414 �����ΪԲ��ʱ �׿����Ϊ�뾶�����²���Ϊ0
            Com_element_main.Ca_cl(i).SideSlopes,&
            Com_element_main.Ca_cl(i).InInvertElev,&
            Com_element_main.Ca_cl(i).OutInvertElev,&
            Com_element_main.Ca_cl(i).InTopElev,&
            Com_element_main.Ca_cl(i).OutTopElev,&
            Com_element_main.Ca_cl(i).Roughness,&
            Com_element_main.Ca_cl(i).CroSecTYPE
    END DO

    !!!///���뽥��β���-Ԫ������7
    IF(Com_element_main.Areachange_NUM/=0) THEN
        READ(12,*)
        DO i=1,Com_element_main.Areachange_NUM
            READ(12,*)Com_element_main.Arch_cl(i).Serial_NUM,&
                Com_element_main.Arch_cl(i).Name,&
                Com_element_main.Arch_cl(i).Incoord,&
                Com_element_main.Arch_cl(i).Outcoord,&
                Com_element_main.Arch_cl(i).Locallosscoeff,&
                Com_element_main.Arch_cl(i).Additionlosscoeff,&
                Com_element_main.Arch_cl(i).SideSlopes
        END DO
    END IF

    !!!///���뵹��������-Ԫ������10
    IF(Com_element_main.Invertedsiphon_NUM/=0)THEN
        READ(12,*)
        DO i=1,Com_element_main.Invertedsiphon_NUM
            READ(12,*)Com_element_main.Insi_cl(i).Serial_NUM,&
                Com_element_main.Insi_cl(i).Name,&
                Com_element_main.Insi_cl(i).Incoord,&
                Com_element_main.Insi_cl(i).Outcoord,&
                Com_element_main.Insi_cl(i).ParallelNum,&
                Com_element_main.Insi_cl(i).BottomWidth,&
                Com_element_main.Insi_cl(i).Height,&
                Com_element_main.Insi_cl(i).Length,&
                Com_element_main.Insi_cl(i).InInvertElev,&
                Com_element_main.Insi_cl(i).OutInvertElev,&
                Com_element_main.Insi_cl(i).Roughness,&
                Com_element_main.Insi_cl(i).Insi_coeff
        END DO
    END IF

    !!!///������������-Ԫ������11
    IF (Com_element_main.Bridge_NUM/=0) THEN
        READ(12,*)
        DO i=1,Com_element_main.Bridge_NUM
            READ(12,*)Com_element_main.Brid_cl(i).Serial_NUM,&
                Com_element_main.Brid_cl(i).Name,&
                Com_element_main.Brid_cl(i).Incoord,&
                Com_element_main.Brid_cl(i).Outcoord,&
                Com_element_main.Brid_cl(i).BottomWidth,&
                Com_element_main.Brid_cl(i).SideSlopes,&
                Com_element_main.Brid_cl(i).InInvertElev,&
                Com_element_main.Brid_cl(i).OutInvertElev,&
                Com_element_main.Brid_cl(i).Roughness,&
                Com_element_main.Brid_cl(i).HinderWidth
        END DO
    END IF

    !!!///�������բ����-Ԫ������6
    IF (Com_element_main.Controlgate_NUM/=0) THEN
        READ(12,*)
        DO i=1,Com_element_main.Controlgate_NUM
            READ(12,*)Com_element_main.Cogate_cl(i).Serial_NUM,&
                Com_element_main.Cogate_cl(i).Name,&
                Com_element_main.Cogate_cl(i).Incoord,&
                Com_element_main.Cogate_cl(i).Outcoord,&
                Com_element_main.Cogate_cl(i).ParallelNum,&
                Com_element_main.Cogate_cl(i).SingleWidth,&
                Com_element_main.Cogate_cl(i).InvertElev,&
                Com_element_main.Cogate_cl(i).RUNNum,&
                Com_element_main.Cogate_cl(i).Discoeff,&
                Com_element_main.Cogate_cl(i).Lenth
        END DO
    END IF

    !!!///�����ˮ�ڲ���-Ԫ������4
    IF (Com_element_main.SideOutFlow_NUM/=0) THEN
        READ(12,*)
        DO i=1,Com_element_main.SideOutFlow_NUM
            READ(12,*)Com_element_main.SideOF_cl(i).Serial_NUM,&
                Com_element_main.SideOF_cl(i).Name,&
                Com_element_main.SideOF_cl(i).Incoord,&
                Com_element_main.SideOF_cl(i).Outcoord
        END DO
    END IF

    !!!///�����վ����-Ԫ������12
    IF(Com_element_main.PumpST_NUM/=0)THEN
        READ(12,*)
        DO i=1,Com_element_main.PumpST_NUM
            READ(12,*)PumpST_cl(i).Serial_NUM,&
                Com_element_main.PumpST_cl(i).Name,&
                Com_element_main.PumpST_cl(i).Incoord,&
                Com_element_main.PumpST_cl(i).Outcoord,&
                Com_element_main.PumpST_cl(i).InvertElev,&
                Com_element_main.PumpST_cl(i).UNIT_NUM
            ALLOCATE(Com_element_main.PumpST_cl(i).Pump_Type(Com_element_main.PumpST_cl(i).UNIT_NUM))
            READ(12,*)(Com_element_main.PumpST_cl(i).Pump_Type(j),j=1,Com_element_main.PumpST_cl(i).UNIT_NUM)
        END DO
    ENDIF

    !**************************
    !******** �豸��Ϣ ********
    !**************************
    !READ(12,*)!================================================================================
    !!!///�豸��Ϣ����ˮ�û���/ˮ�ֻ����顢������Ϣ��
    !!!///���벻ͬ����ˮ�ò���
    IF(Com_element_main.PumpST_NUM/=0)THEN
        READ(12,*) !ˮ�û�������
        READ(12,*)Com_element_main.Pump_type_NUM
        Pump_type_NUM=Com_element_main.Pump_type_NUM
        CALL Com_element_Auto_Update(Com_element_main,"Pump_type_NUM")
        ALLOCATE(Pump_cl(Pump_type_NUM))
        READ(12,*)
        DO i=1,Pump_type_NUM
            READ(12,*)Com_element_main.Pump_cl(i).Serial_NUM,&
                Com_element_main.Pump_cl(i).Name
            READ(12,*)Com_element_main.Pump_cl(i).Coefficient_NUM
            ALLOCATE(Com_element_main.Pump_cl(i).Bladeangle_Serial(Com_element_main.Pump_cl(i).Coefficient_NUM))
            ALLOCATE(Com_element_main.Pump_cl(i).Coefficient_A(Com_element_main.Pump_cl(i).Coefficient_NUM))
            ALLOCATE(Com_element_main.Pump_cl(i).Coefficient_B(Com_element_main.Pump_cl(i).Coefficient_NUM))
            ALLOCATE(Com_element_main.Pump_cl(i).Coefficient_C(Com_element_main.Pump_cl(i).Coefficient_NUM))
            DO j=1,Com_element_main.Pump_cl(i).Coefficient_NUM
                READ(12,*)Com_element_main.Pump_cl(i).Bladeangle_Serial(j),&
                    Com_element_main.Pump_cl(i).Coefficient_A(j),&
                    Com_element_main.Pump_cl(i).Coefficient_B(j),&
                    Com_element_main.Pump_cl(i).Coefficient_C(j)
            END DO
        END DO
    ENDIF
    !!!!!!!Ϊ��վϵ������ռ�
    IF(Com_element_main.Pump_type_NUM/=0.and.Com_element_main.PumpST_NUM/=0)THEN
        DO i=1,Com_element_main.PumpST_NUM
            PS_Co_NUM=1
            DO j=1,Com_element_main.PumpST_cl(i).UNIT_NUM!��վ��ÿ̨�����ת�������
                PS_Co_NUM=PS_Co_NUM*Com_element_main.Pump_cl(Com_element_main.PumpST_cl(i).Pump_Type(j)).Coefficient_NUM
            END DO
            ALLOCATE(Com_element_main.PumpST_cl(i).Coefficient_A(PS_Co_NUM))
            Com_element_main.PumpST_cl(i).Coefficient_A=0.0/0.0
            ALLOCATE(Com_element_main.PumpST_cl(i).Coefficient_B(PS_Co_NUM))
            Com_element_main.PumpST_cl(i).Coefficient_B=0.0/0.0
            ALLOCATE(Com_element_main.PumpST_cl(i).Coefficient_C(PS_Co_NUM))
            Com_element_main.PumpST_cl(i).Coefficient_C=0.0/0.0
        END DO
    ENDIF


    !TU����һ�������ͳ����Ҫ�Ķ������
    !TU����˼·���ǣ�ÿ��n�������Ĳ�������������ô����Ҫn+1������
    !TU�м��һ�����⽨���������ǹ���������Ȼ������m�������Ĳ�������������ô�͵��ټ�m+1�����棻�Դ�����
    irr_num=0
    section_num=0
    irrseries_num=0
    sum_of_irr=0!��������������
    Do i=1,Com_element_main.Element_NUM-1!ͳ���Ƿ�Ϊ�����Ĳ���������
        IF((Com_element_main.Ele_Rel(i,2)==5).and.(Com_element_main.Ca_cl(Com_element_main.Ele_Rel(i,3)).CroSecTYPE==2)) then
            sum_of_irr =sum_of_irr+1
            Com_element_main.Ca_cl(Com_element_main.Ele_Rel(i,3)).irr_Serial_NUM=sum_of_irr
            irr_num= irr_num+1
            IF (.not.((Com_element_main.Ele_Rel(i+1,2)==5).and.(Com_element_main.Ca_cl(Com_element_main.Ele_Rel(i+1,3)).CroSecTYPE==2))) then!TU �������ǲ���������������һ���������ǲ�������������ô����������������+1
                irrseries_num=irrseries_num+1
                section_num = section_num +irr_num + 1
                irr_num_mat(irrseries_num)=irr_num
                irr_num = 0!һ��һ���β��������ͳ����ĩβ���ͽ���������������0
            END IF
        END IF
    END DO
    !����ͳ�����һ���Ƿ�Ϊ����������
    IF(Com_element_main.Ele_Rel(Com_element_main.Element_NUM,2)==5)then
        IF(Com_element_main.Ca_cl(Com_element_main.Ele_Rel(Com_element_main.Element_NUM,3)).CroSecTYPE==2) then!ֻҪ���һ���ǲ�������������ô������Σ��ֲ�irr_num+1,����irrseries_num+1
            !��Ϊ���ǰһ���ǲ�������棬��ô�����Ǹ��ڲ�if��û������irr_num�ͱ��������ˣ�������һ�ε�irrseries_num��û�м�
            !�����ǰһ����������ʲô�������ôirr_numԭ����0����Ҫ�¼�һ��irrseries_num
            sum_of_irr =sum_of_irr+1
            Com_element_main.Ca_cl(Com_element_main.Ele_Rel(i,3)).irr_Serial_NUM=sum_of_irr
            irr_num= irr_num+1
            irrseries_num=irrseries_num+1
            section_num =section_num+irr_num+1
            irr_num_mat(irrseries_num)=irr_num
        END IF
    END IF

    !TU��ȡ�̶߳ԣ���������
    Allocate(SectionName(section_num))
    Allocate(CDS(section_num))
    DO i=1,section_num
        READ(12,*) SectionName(i),CDS(i)  !�����ţ��̶߳�TU,ע����Ҫ���̶߳Եĸ���д�ڶ�������֮��
        DO j=1,CDS(i)
            READ(12,*) GCD_temp(1,j,i),GCD_temp(2,j,i)
        END DO
    END DO

    !TU���̶߳Ը�ֵ��������������ÿ�������������������
    kk=0
    j=1!��ʾĿǰ���ڵڼ���β���������
    do i = 1,Com_element_main.Canal_NUM
        if (Com_element_main.Ca_cl(i).CroSecTYPE==2) then!������в���������
            kk=kk+1
            if (kk-j+1>sum(irr_num_mat(1:j))) then
                kk=kk+1
                j=j+1!������һ��β���������
            end if
            Com_element_main.Ca_cl(i).CDS_in =CDS(kk)
            Com_element_main.Ca_cl(i).CDS_out =CDS(kk+1)
            Com_element_main.Ca_cl(i).GCD_in(1:2,1:CDS(kk)) = GCD_temp(1:2,1:CDS(kk),kk)
            Com_element_main.Ca_cl(i).GCD_out(1:2,1:CDS(kk+1)) = GCD_temp(1:2,1:CDS(kk+1),kk+1)
            Com_element_main.Ca_cl(i).ININVERTELEV=minval(Com_element_main.Ca_cl(i).GCD_in(2,1:CDS(kk)))      !TuֻҪ�ǲ�������棬����ڵĵ׸̡߳����̣߳��׿�ͱ��¶�����д0
            Com_element_main.Ca_cl(i).OUTINVERTELEV=minval(Com_element_main.Ca_cl(i).GCD_out(2,1:CDS(kk+1)))
            Com_element_main.Ca_cl(i).INTOPELEV = maxval(Com_element_main.Ca_cl(i).GCD_in(2,1:CDS(kk)))
            Com_element_main.Ca_cl(i).OUTTOPELEV = maxval(Com_element_main.Ca_cl(i).GCD_out(2,1:CDS(kk+1)))
        end if
    end do

    Deallocate(SectionName)
    Deallocate(CDS)
    CLOSE(12)


    !**************************
    !******** ��ʼ���� ********
    !**************************

    OPEN(UNIT=13,FILE="INPUTFILE/INITIAL_STATUS.txt")
    !���γ�ʼ����,���γ�ʼˮλ
    READ(13,*)
    READ(13,*)!================================================================================
    READ(13,*)
    READ(13,*)Boundary_in.UpQ(1),Boundary_in.DownZ(1)    !���γ�ʼ����,���γ�ʼˮ��



    DO i=1,Boundary_in.NUM !�߽�����
        Boundary_in.DownZ(i)=Boundary_in.DownZ(1)
    END DO

    !DO i=1,Boundary_in.NUM ! �߽�����
    !    IF (i == 1) THEN
    !        Boundary_in.UpQ(i) = Boundary_in.UpQ(1) ! ���ֵ�һ��ʱ�̵�ֵ����
    !    ELSE
    !        Boundary_in.UpQ(i) = MAX(Boundary_in.UpQ(i-1) - inGate_changelimit*Boundary_in.Step/60 , 0.0) ! �ӵڶ�ʱ�̿�ʼÿ�μ���5����ͱ���Ϊ0
    !    END IF
    !END DO
    DO i=1,Boundary_in.NUM !�߽�����
        Boundary_in.UpQ(i)=Boundary_in.UpQ(1)
    END DO



    READ(13,*)!================================================================================

    !բǰ��ʼˮλ������բ��
    IF(Com_element_main.Controlgate_NUM/=0)THEN
        READ(13,*)
        DO i=1,Com_element_main.Controlgate_NUM
            READ(13,*)Com_element_main.Cogate_cl(i).Initial_gate_Zup
        END DO

        READ(13,*)!================================================================================
        READ(13,*)!��ʼբ�ſ���
        READ(13,*) Boundary_in.Opendegree(1,1)
        DO j=1,Boundary_in.NUM
            Boundary_in.Opendegree(1,j)=Boundary_in.Opendegree(1,1)
        END DO
        Boundary_in.Opendegree(:,:)=Boundary_in.Opendegree(:,:)/1000
    END IF
    CLOSE(13)



    !**************************
    !******** ���ع��� ********
    !**************************

    !/////����բ���ع���
    !IF(Boundary_in.Gate_NUM/=0)THEN
    !OPEN(UNIT=15,FILE="INPUTFILE/SluiceGate_BY.txt")
    !�Խ���բ���ȵ���Ҫ���д�����Ҫ�����ܹ�����ÿ��բ�ŵĿ��ȣ��Լ�ʵ�ʹ�ˮ��ȵȡ�
    !READ(15,*)  !բ�ŵ��ر߽�
    !READ(15,*)  !բ�ű��
    !READ(15,*) Boundary_in.Opendegree(1,1)
    !DO j=1,Boundary_in.NUM
    !    Boundary_in.Opendegree(1,j)=Boundary_in.Opendegree(1,1)
    !END DO
    !Boundary_in.Opendegree(:,:)=Boundary_in.Opendegree(:,:)/1000
    !CLOSE(15)
    !END IF

    !/////��վ���ع���
    IF(Boundary_in.PumpST_NUM/=0)THEN
        OPEN(UNIT=16,FILE="INPUTFILE/PumpingStation_BY.txt")
        READ(16,*)  !��վ���ر߽�
        READ(16,*)  !��վ���
        READ(16,*)  !������
        DO i=1,Boundary_in.PumpST_NUM
            ALLOCATE(Boundary_in.Bladeangle_cl(i).Each_PS(Com_element_main.PumpST_cl(i).UNIT_NUM,Boundary_in.NUM))
        END DO
        DO k=1,Boundary_in.NUM
            READ(16,*)((Boundary_in.Bladeangle_cl(i).Each_PS(j,k),j=1,Com_element_main.PumpST_cl(i).UNIT_NUM),i=1,Boundary_in.PumpST_NUM)
        END DO
        CLOSE(16)

    END IF

    !/////��ˮ�ڵ��ع���
    IF(Boundary_in.SideOutFlow_NUM/=0)THEN
        OPEN(UNIT=17,FILE="INPUTFILE/SideOutFlow_BY.txt")
        READ(17,*)  !��ˮ���
        READ(17,*)  !��ˮ�ڱ��
        READ(17,*)(Boundary_in.OutDischarge(i,1),i=1,Boundary_in.SideOutFlow_NUM)

        CLOSE(17)
        DO j=1,Boundary_in.NUM
            DO i=1,Boundary_in.SideOutFlow_NUM
                Boundary_in.OutDischarge(i,j)=Boundary_in.OutDischarge(i,1)
            end DO
        END DO
    ENDIF


    pollut_ConveyDistance=(Com_element_main.Cogate_cl(1).Incoord-Com_element_main.Ca_cl(1).Incoord)*1000-pollut_location





    PumpST_NUM=Com_element_main.PumpST_NUM
    PumpST_cl=Com_element_main.PumpST_cl
    Pump_cl=Com_element_main.Pump_cl


    Bladeangle_cl=Boundary_in.Bladeangle_cl


    END
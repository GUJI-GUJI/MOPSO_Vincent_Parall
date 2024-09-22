    SUBROUTINE Inputfile        !输入文件
    !**************************
    !******** 更新记录 ********
    !**************************
    !优化输入文件框架 by Vincent Zhang 23/5/2022

    USE GLOBAL_VARIABLE

    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_in,Boundary_cal
    !TYPE(Com_element)::Com_element_main
    !TYPE(Config)::Config_main

    REAL::Up_tmp,Down_tmp

    INTEGER::i,j,k,PS_Co_NUM
    INTEGER::irr_num!#tu每一段连续的不规则渠道数统计
    INTEGER::section_num!#TU将要读取的断面数
    INTEGER::irrseries_num!#tu 连续的 不规则渠道 段数，即很多段不规则渠道首尾相接算一段
    INTEGER::irr_num_mat(200)!#tu
    INTEGER,allocatable::CDS(:)!#读取断面时的高程对数
    CHARACTER(len=12),allocatable::SectionName(:)!#tu断面名称
    REAL::GCD_temp(2,3000,200)!#一维二维是高程点，三维是断面序号
    INTEGER::kk,temp!#tu
    INTEGER::sum_of_irr!#tu不规则渠道总数

    !**************************
    !******** 固有变量 ********
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



    !/////配置文件
    OPEN(UNIT=18,FILE="INPUTFILE/Config.txt")
    READ(18,*) !配置
    READ(18,*) !污染物质量
    READ(18,*) pollut_weight
    READ(18,*) !距离上节制闸距离(m)
    READ(18,*) pollut_location
    READ(18,*) !退水闸编号
    READ(18,*) outflow_index
    READ(18,*) !退水闸最大流量
    READ(18,*) outflow_maxQ
    READ(18,*) !退水闸变幅限制
    READ(18,*) outflow_changelimit
    READ(18,*) !上节制闸变幅限制（流量）
    READ(18,*) inGate_changelimit
    READ(18,*) !到达时间权重
    READ(18,*) arrivetime_weight
    READ(18,*) !种群个数
    READ(18,*) populations
    READ(18,*) !迭代次数
    READ(18,*) generations
    CLOSE(18)



    !**************************
    !******** 边界条件 ********
    !**************************
    Boundary_in.Step=300
    Boundary_cal.Step=30
    Boundary_in.Start_Time=0
    Boundary_in.End_Time=24000!暂时手动固定
    Out_DeltaT = Boundary_in.Step
    Boundary_2cal.Step=Boundary_cal.Step
    !Boundary_in.End_Time=Boundary_in.Start_Time+(NP/1)*Boundary_in.Step!这里的2代表了节制闸个数，目前只能手动输入了
    CALL Boundary_Auto_Update(Boundary_in,"NUM")
    !READ(11,*)
    !READ(11,*)Boundary_in.Upsign,Boundary_in.Downsign   !边界条件控制符，0表示流量边界，1表示水位边界，2表示水位流量关系边界（2仅能发生在下游）
    Boundary_in.Upsign=0
    Boundary_in.Downsign=1
    Output_sign=1
    !READ(11,*)
    !READ(11,*)Output_sign            !输出结果控制符，0-只输出节制闸处结果；1-输出各元件断面结果
    !READ(11,*)
    !DO i=1,Boundary_in.NUM !读边界条件
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

    !//////计算元件基本信息（断面、里程、长度）
    OPEN(UNIT=12,FILE="INPUTFILE/COM-ELEMENT.txt")
    READ(12,*) !计算元件基本参数及相互连接关系文件
    READ(12,*)!================================================================================
    READ(12,*) !元件个数
    READ(12,*)Com_element_main.Element_NUM

    READ(12,*) !各类型个数
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
    !///读入渠道参数-元件代码5
    READ(12,*)
    DO i=1,Com_element_main.Canal_NUM
        READ(12,*)Com_element_main.Ca_cl(i).Serial_NUM,&
            Com_element_main.Ca_cl(i).Name,&
            Com_element_main.Ca_cl(i).Incoord,&
            Com_element_main.Ca_cl(i).Outcoord,&
            Com_element_main.Ca_cl(i).BottomWidth, &              !wjb0414 横断面为圆形时 底宽参数为半径，边坡参数为0
            Com_element_main.Ca_cl(i).SideSlopes,&
            Com_element_main.Ca_cl(i).InInvertElev,&
            Com_element_main.Ca_cl(i).OutInvertElev,&
            Com_element_main.Ca_cl(i).InTopElev,&
            Com_element_main.Ca_cl(i).OutTopElev,&
            Com_element_main.Ca_cl(i).Roughness,&
            Com_element_main.Ca_cl(i).CroSecTYPE
    END DO

    !!!///读入渐变段参数-元件代码7
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

    !!!///读入倒虹吸参数-元件代码10
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

    !!!///读入桥梁参数-元件代码11
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

    !!!///读入节制闸参数-元件代码6
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

    !!!///读入分水口参数-元件代码4
    IF (Com_element_main.SideOutFlow_NUM/=0) THEN
        READ(12,*)
        DO i=1,Com_element_main.SideOutFlow_NUM
            READ(12,*)Com_element_main.SideOF_cl(i).Serial_NUM,&
                Com_element_main.SideOF_cl(i).Name,&
                Com_element_main.SideOF_cl(i).Incoord,&
                Com_element_main.SideOF_cl(i).Outcoord
        END DO
    END IF

    !!!///读入泵站参数-元件代码12
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
    !******** 设备信息 ********
    !**************************
    !READ(12,*)!================================================================================
    !!!///设备信息：如水泵机组/水轮机机组、断面信息等
    !!!///读入不同类型水泵参数
    IF(Com_element_main.PumpST_NUM/=0)THEN
        READ(12,*) !水泵机组类型
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
    !!!!!!!为泵站系数分配空间
    IF(Com_element_main.Pump_type_NUM/=0.and.Com_element_main.PumpST_NUM/=0)THEN
        DO i=1,Com_element_main.PumpST_NUM
            PS_Co_NUM=1
            DO j=1,Com_element_main.PumpST_cl(i).UNIT_NUM!泵站内每台机组的转角数相乘
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


    !TU以下一大段用于统计需要的断面个数
    !TU具体思路就是，每有n段连续的不规则渠道，那么就需要n+1个断面
    !TU中间接一个任意建筑物甚至是规则渠道，然后再来m段连续的不规则渠道，那么就得再加m+1个断面；以此类推
    irr_num=0
    section_num=0
    irrseries_num=0
    sum_of_irr=0!不规则渠道总数
    Do i=1,Com_element_main.Element_NUM-1!统计是否为连续的不规则渠道
        IF((Com_element_main.Ele_Rel(i,2)==5).and.(Com_element_main.Ca_cl(Com_element_main.Ele_Rel(i,3)).CroSecTYPE==2)) then
            sum_of_irr =sum_of_irr+1
            Com_element_main.Ca_cl(Com_element_main.Ele_Rel(i,3)).irr_Serial_NUM=sum_of_irr
            irr_num= irr_num+1
            IF (.not.((Com_element_main.Ele_Rel(i+1,2)==5).and.(Com_element_main.Ca_cl(Com_element_main.Ele_Rel(i+1,3)).CroSecTYPE==2))) then!TU 本渠道是不规则渠道，且下一个渠道不是不规则渠道，那么连续不规则渠道段+1
                irrseries_num=irrseries_num+1
                section_num = section_num +irr_num + 1
                irr_num_mat(irrseries_num)=irr_num
                irr_num = 0!一旦一整段不规则断面统计至末尾，就将不规则渠道数归0
            END IF
        END IF
    END DO
    !单独统计最后一段是否为不规则渠道
    IF(Com_element_main.Ele_Rel(Com_element_main.Element_NUM,2)==5)then
        IF(Com_element_main.Ca_cl(Com_element_main.Ele_Rel(Com_element_main.Element_NUM,3)).CroSecTYPE==2) then!只要最后一段是不规则渠道，那么无论如何，局部irr_num+1,总体irrseries_num+1
            !因为如果前一段是不规则断面，那么上面那个内部if就没触发，irr_num就保留下来了，并且上一段的irrseries_num还没有加
            !而如果前一段是其他的什么建筑物，那么irr_num原本是0，需要新加一段irrseries_num
            sum_of_irr =sum_of_irr+1
            Com_element_main.Ca_cl(Com_element_main.Ele_Rel(i,3)).irr_Serial_NUM=sum_of_irr
            irr_num= irr_num+1
            irrseries_num=irrseries_num+1
            section_num =section_num+irr_num+1
            irr_num_mat(irrseries_num)=irr_num
        END IF
    END IF

    !TU读取高程对，储存起来
    Allocate(SectionName(section_num))
    Allocate(CDS(section_num))
    DO i=1,section_num
        READ(12,*) SectionName(i),CDS(i)  !断面编号，高程对TU,注：需要将高程对的个数写在断面名字之后
        DO j=1,CDS(i)
            READ(12,*) GCD_temp(1,j,i),GCD_temp(2,j,i)
        END DO
    END DO

    !TU将高程对赋值给不规则渠道，每个渠道由两个断面控制
    kk=0
    j=1!表示目前正在第几大段不规则渠道
    do i = 1,Com_element_main.Canal_NUM
        if (Com_element_main.Ca_cl(i).CroSecTYPE==2) then!针对所有不规则渠道
            kk=kk+1
            if (kk-j+1>sum(irr_num_mat(1:j))) then
                kk=kk+1
                j=j+1!进入下一大段不规则渠道
            end if
            Com_element_main.Ca_cl(i).CDS_in =CDS(kk)
            Com_element_main.Ca_cl(i).CDS_out =CDS(kk+1)
            Com_element_main.Ca_cl(i).GCD_in(1:2,1:CDS(kk)) = GCD_temp(1:2,1:CDS(kk),kk)
            Com_element_main.Ca_cl(i).GCD_out(1:2,1:CDS(kk+1)) = GCD_temp(1:2,1:CDS(kk+1),kk+1)
            Com_element_main.Ca_cl(i).ININVERTELEV=minval(Com_element_main.Ca_cl(i).GCD_in(2,1:CDS(kk)))      !Tu只要是不规则断面，出入口的底高程、顶高程，底宽和边坡都可以写0
            Com_element_main.Ca_cl(i).OUTINVERTELEV=minval(Com_element_main.Ca_cl(i).GCD_out(2,1:CDS(kk+1)))
            Com_element_main.Ca_cl(i).INTOPELEV = maxval(Com_element_main.Ca_cl(i).GCD_in(2,1:CDS(kk)))
            Com_element_main.Ca_cl(i).OUTTOPELEV = maxval(Com_element_main.Ca_cl(i).GCD_out(2,1:CDS(kk+1)))
        end if
    end do

    Deallocate(SectionName)
    Deallocate(CDS)
    CLOSE(12)


    !**************************
    !******** 初始条件 ********
    !**************************

    OPEN(UNIT=13,FILE="INPUTFILE/INITIAL_STATUS.txt")
    !上游初始流量,下游初始水位
    READ(13,*)
    READ(13,*)!================================================================================
    READ(13,*)
    READ(13,*)Boundary_in.UpQ(1),Boundary_in.DownZ(1)    !上游初始流量,下游初始水深



    DO i=1,Boundary_in.NUM !边界条件
        Boundary_in.DownZ(i)=Boundary_in.DownZ(1)
    END DO

    !DO i=1,Boundary_in.NUM ! 边界条件
    !    IF (i == 1) THEN
    !        Boundary_in.UpQ(i) = Boundary_in.UpQ(1) ! 保持第一个时刻的值不变
    !    ELSE
    !        Boundary_in.UpQ(i) = MAX(Boundary_in.UpQ(i-1) - inGate_changelimit*Boundary_in.Step/60 , 0.0) ! 从第二时刻开始每次减少5，最低保持为0
    !    END IF
    !END DO
    DO i=1,Boundary_in.NUM !边界条件
        Boundary_in.UpQ(i)=Boundary_in.UpQ(1)
    END DO



    READ(13,*)!================================================================================

    !闸前初始水位（节制闸）
    IF(Com_element_main.Controlgate_NUM/=0)THEN
        READ(13,*)
        DO i=1,Com_element_main.Controlgate_NUM
            READ(13,*)Com_element_main.Cogate_cl(i).Initial_gate_Zup
        END DO

        READ(13,*)!================================================================================
        READ(13,*)!初始闸门开度
        READ(13,*) Boundary_in.Opendegree(1,1)
        DO j=1,Boundary_in.NUM
            Boundary_in.Opendegree(1,j)=Boundary_in.Opendegree(1,1)
        END DO
        Boundary_in.Opendegree(:,:)=Boundary_in.Opendegree(:,:)/1000
    END IF
    CLOSE(13)



    !**************************
    !******** 调控过程 ********
    !**************************

    !/////节制闸调控过程
    !IF(Boundary_in.Gate_NUM/=0)THEN
    !OPEN(UNIT=15,FILE="INPUTFILE/SluiceGate_BY.txt")
    !对节制闸开度等需要进行处理，主要包括能够设置每个闸门的开度，以及实际过水宽度等。
    !READ(15,*)  !闸门调控边界
    !READ(15,*)  !闸门编号
    !READ(15,*) Boundary_in.Opendegree(1,1)
    !DO j=1,Boundary_in.NUM
    !    Boundary_in.Opendegree(1,j)=Boundary_in.Opendegree(1,1)
    !END DO
    !Boundary_in.Opendegree(:,:)=Boundary_in.Opendegree(:,:)/1000
    !CLOSE(15)
    !END IF

    !/////泵站调控过程
    IF(Boundary_in.PumpST_NUM/=0)THEN
        OPEN(UNIT=16,FILE="INPUTFILE/PumpingStation_BY.txt")
        READ(16,*)  !泵站调控边界
        READ(16,*)  !泵站编号
        READ(16,*)  !机组编号
        DO i=1,Boundary_in.PumpST_NUM
            ALLOCATE(Boundary_in.Bladeangle_cl(i).Each_PS(Com_element_main.PumpST_cl(i).UNIT_NUM,Boundary_in.NUM))
        END DO
        DO k=1,Boundary_in.NUM
            READ(16,*)((Boundary_in.Bladeangle_cl(i).Each_PS(j,k),j=1,Com_element_main.PumpST_cl(i).UNIT_NUM),i=1,Boundary_in.PumpST_NUM)
        END DO
        CLOSE(16)

    END IF

    !/////分水口调控过程
    IF(Boundary_in.SideOutFlow_NUM/=0)THEN
        OPEN(UNIT=17,FILE="INPUTFILE/SideOutFlow_BY.txt")
        READ(17,*)  !分水情况
        READ(17,*)  !分水口编号
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
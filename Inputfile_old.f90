SUBROUTINE Inputfile
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    INTEGER::i,j,k,PS_Co_NUM
    
    !**************************
    !******** 固有变量 ********
    !**************************
    Theta=0.6
    Phi=0.5
    Psi=0.6
    PsiR=0.5
    !RoughnessCoeff=0.02
    NUMMember=1
    
    
    !**************************
    !******** 边界条件 ********
    !**************************
    
    OPEN(UNIT=11,FILE="INPUTFILE/BOUNDARY.txt")
        READ(11,*)
        READ(11,*)
        READ(11,*)Start_Time,End_Time,Boundary_Step,DeltaT,Controltime!输出时间步长 
        Boundary_NUM=(End_Time-Start_Time)/Boundary_Step+1
        READ(11,*)
        READ(11,*)Boundary_sign_UP,Boundary_sign_DOWN   !边界条件控制符，0表示流量边界，1表示水位边界，2表示水位流量关系边界（2仅能发生在下游）
        ALLOCATE(R_Upstream(Boundary_NUM))
        ALLOCATE(R_Downstream(Boundary_NUM))
        READ(11,*)
        DO i=1,Boundary_NUM !读边界条件
            READ(11,*)R_Upstream(i),R_Downstream(i)
        END DO
    CLOSE(11)
    !边界条件的插值计算
    IF(Boundary_Step/=DeltaT)THEN
        Calculate_NUM=(End_Time-Start_Time)/DeltaT+1
        ALLOCATE(Upstream(Calculate_NUM))
	    ALLOCATE(Downstream(Calculate_NUM))
        k=Boundary_Step/DeltaT
        DO j=1,Calculate_NUM-2
            Upstream(j)=R_Upstream(j/k+1)+(MOD(j,k)-1)*(R_Upstream(j/k+2)-R_Upstream(j/k+1))/k
            Downstream(j)=R_Downstream(j/k+1)+(MOD(j,k)-1)*(R_Downstream(j/k+2)-R_Downstream(j/k+1))/k
        END DO
        Upstream(Calculate_NUM-1)=R_Upstream(Boundary_NUM)-(R_Upstream(Boundary_NUM)-R_Upstream(Boundary_NUM-1))/k
        Downstream(Calculate_NUM-1)=R_Downstream(Boundary_NUM)-(R_Downstream(Boundary_NUM)-R_Downstream(Boundary_NUM-1))/k
        Upstream(Calculate_NUM)=R_Upstream(Boundary_NUM)
        Downstream(Calculate_NUM)=R_Downstream(Boundary_NUM)
    ELSEIF(Boundary_Step==DeltaT)THEN
         Calculate_NUM=(End_Time-Start_Time)/DeltaT+1
        ALLOCATE(Upstream(Calculate_NUM))
	    ALLOCATE(Downstream(Calculate_NUM))
        k=Boundary_Step/DeltaT
        DO j=1,Calculate_NUM
            Upstream(j)=R_Upstream(j)
            Downstream(j)=R_Downstream(j)
        END DO
    END IF
    
    
    !**************************
    !******** 计算元件 ********
    !**************************
    
    !//////计算元件基本信息（断面、里程、长度）
    OPEN(UNIT=12,FILE="INPUTFILE/COM-ELEMENT.txt")
        READ(12,*) !计算元件基本参数及相互连接关系文件
        READ(12,*)!================================================================================
        READ(12,*) !元件个数
        READ(12,*)Element_NUM
        READ(12,*) !各类型个数
        READ(12,*)
        READ(12,*)Canal_NUM,Areachange_NUM,Invertedsiphon_NUM,Controlgate_NUM,SideOutFlow_NUM,Bridge_NUM,PumpST_NUM
        READ(12,*)!================================================================================
        READ(12,*)
        ALLOCATE(Ele_Rel(Element_NUM,3))
        DO i=1,Element_NUM
            READ(12,*)(Ele_Rel(i,j),j=1,3)
        END DO
        READ(12,*)!================================================================================
        !///读入渠道参数-元件代码5
        ALLOCATE(Ca_cl(Canal_NUM))
        READ(12,*)
        DO i=1,Canal_NUM
            READ(12,*)Ca_cl(i)%Serial_NUM,Ca_cl(i)%Name,Ca_cl(i)%Incoord,Ca_cl(i)%Outcoord,Ca_cl(i)%BottomWidth, &              !wjb0414 横断面为圆形时 底宽参数为半径，边坡参数为0
            Ca_cl(i)%SideSlopes,Ca_cl(i)%InInvertElev,Ca_cl(i)%OutInvertElev,Ca_cl(i)%InTopElev,Ca_cl(i)%OutTopElev,Ca_cl(i)%Roughness, Ca_cl(i)%CroSecTYPE 
        END DO
        
        !!!///读入渐变段参数-元件代码7
        ALLOCATE(Arch_cl(Areachange_NUM))
        IF (Areachange_NUM/=0) THEN
            READ(12,*)
            DO i=1,Areachange_NUM
                 READ(12,*)Arch_cl(i)%Serial_NUM,Arch_cl(i)%Name,Arch_cl(i)%Incoord,Arch_cl(i)%Outcoord,Arch_cl(i)%Locallosscoeff,Arch_cl(i)%Additionlosscoeff,Arch_cl(i)%SideSlopes
            END DO
        END IF
        
        !!!///读入倒虹吸参数-元件代码10
        ALLOCATE(Insi_cl(Invertedsiphon_NUM))
        IF (Invertedsiphon_NUM/=0) THEN
            READ(12,*)
            DO i=1,Invertedsiphon_NUM
                READ(12,*)Insi_cl(i)%Serial_NUM,Insi_cl(i)%Name,Insi_cl(i)%Incoord,Insi_cl(i)%Outcoord,Insi_cl(i)%ParallelNum,Insi_cl(i)%BottomWidth,Insi_cl(i)%Height,Insi_cl(i)%Length,Insi_cl(i)%InInvertElev,Insi_cl(i)%OutInvertElev,Insi_cl(i)%Roughness,Insi_cl(i)%Insi_coeff
            END DO 
        END IF
        
        !!!///读入桥梁参数-元件代码11
        ALLOCATE(Brid_cl(Bridge_NUM))
        IF (Bridge_NUM/=0) THEN
            READ(12,*)
            DO i=1,Bridge_NUM
                READ(12,*)Brid_cl(i)%Serial_NUM,Brid_cl(i)%Name,Brid_cl(i)%Incoord,Brid_cl(i)%Outcoord,Brid_cl(i)%BottomWidth,Brid_cl(i)%SideSlopes,Brid_cl(i)%InInvertElev,Brid_cl(i)%OutInvertElev,Brid_cl(i)%Roughness,Brid_cl(i)%HinderWidth
            END DO
        END IF
        
        !!!///读入节制闸参数-元件代码6
        ALLOCATE(Cogate_cl(Controlgate_NUM))
        IF (Controlgate_NUM/=0) THEN
            READ(12,*)
            DO i=1,Controlgate_NUM
                !READ(12,*)Cogate_cl(i)%Serial_NUM,Cogate_cl(i)%Name,Cogate_cl(i)%Incoord,Cogate_cl(i)%Outcoord,Cogate_cl(i)%ParallelNum,Cogate_cl(i)%SingleWidth,Cogate_cl(i)%InvertElev,Cogate_cl(i)%RUNNum,Cogate_cl(i)%G_A,Cogate_cl(i)%G_B,Cogate_cl(i)%G_C
                READ(12,*)Cogate_cl(i)%Serial_NUM,Cogate_cl(i)%Name,Cogate_cl(i)%Incoord,Cogate_cl(i)%Outcoord,Cogate_cl(i)%ParallelNum,Cogate_cl(i)%SingleWidth,Cogate_cl(i)%InvertElev,Cogate_cl(i)%RUNNum,Cogate_cl(i)%Discoeff,Cogate_cl(i)%Lenth
            END DO
        END IF
        
        !!!///读入分水口参数-元件代码4
        ALLOCATE(SideOF_cl(SideOutFlow_NUM))
        IF (SideOutFlow_NUM/=0) THEN
            READ(12,*)
            DO i=1,SideOutFlow_NUM
                READ(12,*)SideOF_cl(i)%Serial_NUM,SideOF_cl(i)%Name,SideOF_cl(i)%Incoord,SideOF_cl(i)%Outcoord
            END DO
        END IF
        
        !!!///读入泵站参数-元件代码12
        ALLOCATE(PumpST_cl(PumpST_NUM))
        IF(PumpST_NUM/=0)THEN
            READ(12,*)
            DO i=1,PumpST_NUM
                READ(12,*)PumpST_cl(i)%Serial_NUM,PumpST_cl(i)%Name,PumpST_cl(i)%Incoord,PumpST_cl(i)%Outcoord,PumpST_cl(i)%InvertElev,PumpST_cl(i)%UNIT_NUM
                ALLOCATE(PumpST_cl(i)%Pump_Type(PumpST_cl(i)%UNIT_NUM))
                READ(12,*)(PumpST_cl(i)%Pump_Type(j),j=1,PumpST_cl(i)%UNIT_NUM)
            END DO
        ENDIF
        
    
    !**************************
    !******** 设备信息 ********
    !**************************
        READ(12,*)!================================================================================
        !!!///设备信息：如水泵机组/水轮机机组、断面信息等
        !!!///读入不同类型水泵参数
        IF(PumpST_NUM/=0)THEN
            READ(12,*) !水泵机组类型
            READ(12,*)Pump_type_NUM
            ALLOCATE(Pump_cl(Pump_type_NUM))
            READ(12,*)
            DO i=1,Pump_type_NUM
                READ(12,*)Pump_cl(i)%Serial_NUM,Pump_cl(i)%Name
                READ(12,*)Pump_cl(i)%Coefficient_NUM
                ALLOCATE(Pump_cl(i)%Bladeangle_Serial(Pump_cl(i)%Coefficient_NUM))
                ALLOCATE(Pump_cl(i)%Coefficient_A(Pump_cl(i)%Coefficient_NUM))
                ALLOCATE(Pump_cl(i)%Coefficient_B(Pump_cl(i)%Coefficient_NUM))
                ALLOCATE(Pump_cl(i)%Coefficient_C(Pump_cl(i)%Coefficient_NUM))
                DO j=1,Pump_cl(i)%Coefficient_NUM
                    READ(12,*)Pump_cl(i)%Bladeangle_Serial(j),Pump_cl(i)%Coefficient_A(j),Pump_cl(i)%Coefficient_B(j),Pump_cl(i)%Coefficient_C(j)
                END DO
            END DO
        ENDIF
        !!!!!!!为泵站系数分配空间
        IF(Pump_type_NUM/=0.and.PumpST_NUM/=0)THEN
            DO i=1,PumpST_NUM
                PS_Co_NUM=1
                DO j=1,PumpST_cl(i)%UNIT_NUM!泵站内每台机组的转角数相乘
                    PS_Co_NUM=PS_Co_NUM*Pump_cl(PumpST_cl(i)%Pump_Type(j))%Coefficient_NUM
                END DO
                ALLOCATE(PumpST_cl(i)%Coefficient_A(PS_Co_NUM))
                PumpST_cl(i)%Coefficient_A=0.0/0.0
                ALLOCATE(PumpST_cl(i)%Coefficient_B(PS_Co_NUM))
                PumpST_cl(i)%Coefficient_B=0.0/0.0
                ALLOCATE(PumpST_cl(i)%Coefficient_C(PS_Co_NUM))
                PumpST_cl(i)%Coefficient_C=0.0/0.0
            END DO
        ENDIF
        
        READ(12,*)!================================================================================
        !///读取不规则断面测点数据
        READ(12,*)!起点距、高程
        DO i=1,Canal_NUM
            IF(Ca_cl(i)%CroSecTYPE==2) THEN
                READ(12,*)  !断面编号
                Ca_cl(i)%CDS=INT(Ca_cl(i)%BottomWidth)              !wjb0414 横断面为不规则断面时 底宽参数为断面测点数，边坡参数为0
                DO j=1,Ca_cl(i)%CDS
                    READ(12,*) Ca_cl(i)%GCD(1,j),Ca_cl(i)%GCD(2,j) 
                END DO
            END IF
        END DO
        !***********************************************************************************************************************************
    CLOSE(12)
    
    
    !**************************
    !******** 初始条件 ********
    !**************************
    
    OPEN(UNIT=13,FILE="INPUTFILE/INITIAL_STATUS.txt")
        !上游初始流量,下游初始水深
        READ(13,*)
        READ(13,*)!================================================================================
        READ(13,*)
        READ(13,*)Initial_Discharge,D_elevation    !上游初始流量,下游初始水深
        READ(13,*)!================================================================================
        
        !闸前初始水位（节制闸）
        IF(Controlgate_NUM/=0)THEN
            READ(13,*)   
            DO i=1,Controlgate_NUM
                READ(13,*)Cogate_cl(i)%Initial_gate_Zup
            END DO
        END IF
        READ(13,*)!================================================================================
        
        !站前初始水位（泵站）
        IF(PumpST_NUM/=0)THEN
            ALLOCATE(Ini_PumpST_cl(PumpST_NUM))
            READ(13,*)   
            DO i=1,PumpST_NUM
                READ(13,*)Ini_PumpST_cl(i)%Initial_pump_Zup
            END DO
        END IF
    CLOSE(13)
   
    !**************************
    !******** 调控过程 ********
    !**************************
    
    !/////节制闸调控过程
    OPEN(UNIT=21,FILE="INPUTFILE/INNERBOUNDARY.txt")
        !对节制闸开度等需要进行处理，主要包括能够设置每个闸门的开度，以及实际过水宽度等。
        READ(21,*)  !闸门调控边界
        READ(21,*)  !闸门编号
        ALLOCATE(Opendegree_INIT(Controlgate_NUM,1))       
        !///////节制闸初始开度,其余程序生成的开度放在后面的水动力计算中
        DO j=1,1
            READ(21,*)(Opendegree_INIT(i,j),i=1,Controlgate_NUM)
        END DO
    CLOSE(21)

    !/////分水口调控过程
    IF(SideOutFlow_NUM/=0)THEN
        OPEN(UNIT=17,FILE="INPUTFILE/SideOutFlow.txt")
            READ(17,*)  !分水情况
            READ(17,*)  !分水口编号
            ALLOCATE(OutDischarge(SideOutFlow_NUM,Boundary_NUM))
            DO j=1,Boundary_NUM
                READ(17,*)(OutDischarge(i,j),i=1,SideOutFlow_NUM)
            END DO
        CLOSE(17)
        !调控过程加密处理
        IF(Boundary_Step/=DeltaT)THEN
            ALLOCATE(C_OutDischarge(SideOutFlow_NUM,Calculate_NUM))
            k=Boundary_Step/DeltaT
            DO j=1,Calculate_NUM-2
                C_OutDischarge(:,j)=OutDischarge(:,j/k+1)+(MOD(j,k)-1)*(OutDischarge(:,j/k+2)-OutDischarge(:,j/k+1))/k
            END DO
            C_OutDischarge(:,Calculate_NUM-1)=OutDischarge(:,Boundary_NUM)-(OutDischarge(:,Boundary_NUM)-OutDischarge(:,Boundary_NUM-1))/k
            C_OutDischarge(:,Calculate_NUM)=OutDischarge(:,Boundary_NUM)
        ELSEIF(Boundary_Step==DeltaT)THEN
            ALLOCATE(C_OutDischarge(SideOutFlow_NUM,Calculate_NUM))
            k=Boundary_Step/DeltaT
            DO j=1,Calculate_NUM
                C_OutDischarge(:,j)=OutDischarge(:,j)
            END DO
        END IF
    ENDIF
    
    OPEN(UNIT=14,FILE="INPUTFILE/TARGET.txt")   !节制闸目标水位
        READ(14,*)
        ALLOCATE(ZQSW_TARGET(Controlgate_NUM))
        DO i=1,Controlgate_NUM
            READ(14,*)ZQSW_TARGET(i)
!            PRINT *,ZQSW_TARGET(i)
        ENDDO
        READ(14,*)D_WL_TARGET
        READ(14,*)
        ALLOCATE(ST_CON_UP(Controlgate_NUM))
        ALLOCATE(ST_CON_DOWN(Controlgate_NUM))
        DO i=1,Controlgate_NUM
            READ(14,*)ST_CON_UP(i),ST_CON_DOWN(i)    
        ENDDO
        READ(14,*)DS_CON_UP,DS_CON_DOWN !下游水位边界，上限、下限
    CLOSE(14)
    ALLOCATE(BestOpendegree(Controlgate_NUM,Boundary_NUM))
    ALLOCATE(BestST_STATE((End_Time-Start_Time)/3600+1,2*Controlgate_NUM))
    ALLOCATE(BestQControlGate((End_Time-Start_Time)/3600+1,2*Controlgate_NUM))
    ALLOCATE(BestVolume(Controlgate_NUM+1,Calculate_NUM))
END
   

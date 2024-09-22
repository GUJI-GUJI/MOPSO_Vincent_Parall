    SUBROUTINE Dynamic_Cal(Boundary_temp,pollu_Lenth,arrive_time)  !非恒定流计算
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
    INTEGER,ALLOCATABLE::DM_POSTION(:)  !节制闸+首尾断面位置
    INTEGER,ALLOCATABLE::SG_POSTION(:)  !节制闸位置
    REAL::m_coeff    !末端闸门过闸流量系数计算
    REAL::PPQ, VVQ  !水位~流量关系断面系数
    INTEGER::Canal_NUM
    INTEGER::c_element    !计算调整空间步长后参与计算元件的个数

    INTEGER::i,j,k,mm,kk,LLL
    real::c11,c22,c21,hh
    REAL::zzaverage,QQaverage
    INTEGER::si_n   !优化存储时的循环变量
    REAL::random    !随机数
    DOUBLE PRECISION r0
    Integer::pos_sg,DM_sg !记录节制闸位置
    INTEGER::ele_sg    !记录各元件的位置（除了桥墩）
    !稳态各断面水深数组（新）
    REAL,ALLOCATABLE::h_new(:)

    !/////////////非恒定流计算
    !未知量存储
    REAL,ALLOCATABLE::ZZ(:,:)
    REAL,ALLOCATABLE::QQ(:,:)
    REAL,ALLOCATABLE::Pool_Volume(:,:)

    !上一步的一直量存储
    REAL,ALLOCATABLE::ZZOLD(:,:)
    REAL,ALLOCATABLE::QQOLD(:,:)
    REAL,ALLOCATABLE::ELE_Volume(:,:)
    !追赶系数
    REAL,ALLOCATABLE::SS(:)
    REAL,ALLOCATABLE::TT(:)



    Real::GateBottom,GateWidth,BeginBottom
    real::pollu_Lenth
    real::arrive_time
    real::t       !时间变量,单位：min
    real::D  !污染源距分水口500米处距离
    real::S  !计算得到的去染污前端距污染源距离
    real::differ    !误差
    real::D1,D2,w    !迭代计算得出的D
    REAL::v,m,B,h,JJ,g   !系数
    real,allocatable::v_eachTime(:)!每次计算的速度
    real,allocatable::h_eachTime(:)!每次计算的h

    allocate(v_eachTime(Boundary_cal.NUM))
    allocate(h_eachTime(Boundary_cal.NUM))
    !allocate(Opendegree(1,Boundary_cal.NUM))

    v=1.0
    m=1.0    !自由倍数取1
    B=26.5
    !h=6.25
    JJ=0.00001
    g=pollut_weight !污染物质量。该例子为：10t苯酚水，计算单位:g,后续可手动更换污染物质量
    D=pollut_ConveyDistance    !单位:m

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
    !未知量存储
    ALLOCATE(ZZ(c_element+1,NUMMember))
    ALLOCATE(QQ(c_element+1,NUMMember))
    ALLOCATE(ZZOLD(c_element+1,NUMMember))
    ALLOCATE(QQOLD(c_element+1,NUMMember))
    ALLOCATE(SG_POSTION(Com_element_cal.Controlgate_NUM))
    ALLOCATE(DM_POSTION(Com_element_cal.Controlgate_NUM+2))

    ALLOCATE(Dynamic_Result_cal.Pool_Volume(Com_element_cal.Controlgate_NUM+1,Boundary_cal.NUM))
    ALLOCATE(ElE_Volume(c_element,NUMMember))
    ALLOCATE(Dynamic_Result_cal.Total_Volume(Boundary_cal.NUM))
    !///////确定闸门所在的位置//////
    pos_sg=1
    DO j=1,c_element
        IF(Com_element_cal.Ele_Rel(j,2)==6)THEN
            SG_POSTION(pos_sg)=j
            pos_sg=pos_sg+1
        ENDIF
    END DO

    !///////确定闸门+首尾断面所在的位置//////
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
            IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////渠道元件
                ZZ(j,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j,3)).InInvertElev+h_new(j)
                ZZ(j+1,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j,3)).OutInvertElev+h_new(j+1)
            ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
                ZZ(j,mm)=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j,3)).InInvertElev+h_new(j)
                ZZ(j+1,mm)=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j,3)).OutInvertElev+h_new(j+1)
            ENDIF
        END DO
    END DO
    IF(Boundary_cal.Upsign==1)THEN
        ZZ(1,1)=Boundary_cal.UpZ(1)
    ENDIF

    !////////总的计算次数为：Boundary_cal.NUM,在输入文件中计算过了
    si_n=1
    DO k=1,Boundary_cal.NUM
        DO mm=1,NUMMember
            DO i=1,c_element,1
                QQOLD(i,mm)=QQ(i,mm)
                if((Com_element_cal.Ele_Rel(i,2)==5).and.( ZZ(i,mm)<Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(i,3)).InInvertElev+0.1)) then
                    ZZOLD(i,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(i,3)).InInvertElev+0.1
                    write(*,*)"==============warning（追赶法计算结果）==============="
                    write(*,*) "第",k,"次计算的时候露底了,将水深赋值为0.1继续计算"
                    write(*,*) "如果连续出现这里的露底问题，会导致结果崩溃。具体原因暂未可知"
                    write(*,*) "可以尝试将计算步长缩短，或者调试一下糙率"
                    write(*,*)"================================================"
                else
                    ZZOLD(i,mm)=ZZ(i,mm)
                end if
            END DO
            !最后一个断面单独处理
            QQOLD(c_element+1,mm)=QQ(c_element+1,mm)
            if((Com_element_cal.Ele_Rel(c_element,2)==5).and.( ZZ(c_element,mm)<Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(c_element,3)).OutInvertElev+0.1)) then
                ZZOLD(c_element+1,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(c_element,3)).OutInvertElev+0.1
            else
                ZZOLD(c_element+1,mm)=ZZ(c_element+1,mm)
            end if
            !#TU-------------------------修改到此结束
            !追赶系数
            ALLOCATE(SS(c_element+1))
            ALLOCATE(TT(c_element+1))
            ALLOCATE(Condition_cal.PP(c_element+1))
            ALLOCATE(Condition_cal.VV(c_element+1))


            !***************生成系数************
            IF(Boundary_cal.Upsign==0)THEN        !上游边界条件分别选择水位和流量.0表示选择流量边界，1表示上游为水位边界
                !/////////上游边界系数
                Condition_cal.PP(1)=Boundary_cal.UpQ(k)
                QQ(1,mm)=Boundary_cal.UpQ(k)
                Condition_cal.VV(1)=0.0
                !/////////中间元件系数
                DO j=1,c_element
                    !/////先判断类型再调用对应的函数模块，生成相应的系数。
                    IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////渠道元件
                        Call CanalDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==4)THEN     !////分水口元件
                        Call SideOutFlowDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_cal.OutDischarge(Com_element_cal.Ele_Rel(j,3),k),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==6)THEN     !////节制闸元件
                        Call ControlGateDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_temp.Opendegree(1,k),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                        firstgate_Index = j
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==7)THEN  !/////渐变段元件
                        Call AreachangeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
                        Call InvertedsiphonDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==11)THEN !/////桥墩元件
                        Call BridgeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==12)THEN    !/////泵站元件
                        Call PumpingstationDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),C_Bladeangle_cl(Com_element_cal.Ele_Rel(j,3)).Each_PS(:,k),Com_element_cal.PumpST_cl(Com_element_cal.Ele_Rel(j,3)).UNIT_NUM,QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ENDIF
                END DO
                !///////////下游边界条件
                IF(Boundary_cal.Downsign==1)THEN   !///////////水位边界1
                    ZZ(c_element+1,mm)=Boundary_cal.DownZ(k)
                ELSEIF(Boundary_cal.Downsign==2)THEN   !///////////水位~流量关系2    此处的下游边界均使用的DownZ，可能会有问题
                    PPQ=m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))-1.0/(2.0*m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))
                    VVQ=-1.0/(2.0*m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev)))
                    ZZ(c_element+1,mm)=(Condition_cal.PP(c_element+1)-PPQ)/(Condition_cal.VV(c_element+1)-VVQ)
                ENDIF

                !*******回代求解水动力过程********！
                DO j=c_element,1,-1
                    ZZ(j,mm)=SS(j+1)-TT(j+1)*ZZ(j+1,mm)
                    QQ(j+1,mm)=Condition_cal.PP(j+1)-Condition_cal.VV(j+1)*ZZ(j+1,mm)
                END DO
            ELSE    !上游为水位边界Boundary_cal.Upsign==1
                !/////////上游边界系数
                Condition_cal.PP(1)=Boundary_cal.UpZ(k)
                ZZ(1,mm)=Boundary_cal.UpZ(k)
                Condition_cal.VV(1)=0.0
                !/////////中间元件系数
                DO j=1,c_element
                    !/////先判断类型再调用对应的函数模块，生成相应的系数。
                    IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////渠道元件
                        Call CanalDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==4)THEN     !////分水口元件
                        Call SideOutFlowDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_cal.OutDischarge(Com_element_cal.Ele_Rel(j,3),k),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==6)THEN     !////节制闸元件
                        Call ControlGateDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_temp.Opendegree(1,k),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                        firstgate_Index = j
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==7)THEN  !/////渐变段元件
                        Call AreachangeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
                        Call InvertedsiphonDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==11)THEN !/////桥墩元件
                        Call BridgeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==12)THEN    !/////泵站元件
                        Call PumpingstationDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),C_Bladeangle_cl(Com_element_cal.Ele_Rel(j,3)).Each_PS(:,k),Com_element_cal.PumpST_cl(Com_element_cal.Ele_Rel(j,3)).UNIT_NUM,QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ENDIF
                END DO
                !///////////下游边界条件
                IF(Boundary_cal.Downsign==0)THEN   !///////////流量边界0
                    QQ(c_element+1,mm)=Boundary_cal.DownQ(k)
                ELSEIF(Boundary_cal.Downsign==1)THEN   !///////////水位边界1
                    QQ(c_element+1,mm)=(Condition_cal.PP(c_element+1)-Boundary_cal.DownZ(k))/Condition_cal.VV(c_element+1)
                ELSEIF(Boundary_cal.Downsign==2)THEN   !///////////水位~流量关系边界2
                    PPQ=m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))-1.0/(2.0*m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))
                    VVQ=-1.0/(2.0*m_coeff*5*Boundary_cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev)))
                    QQ(c_element+1,mm)=(PPQ-VVQ*Condition_cal.PP(c_element+1))/(1-VVQ*Condition_cal.VV(c_element+1))
                ENDIF
                !*******回代求解水动力过程********！
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

            !***各渠池求总和*******
            DO i=1,Com_element_cal.Controlgate_NUM+1
                Dynamic_Result_cal.Pool_Volume(i,k)=0.0
                DO j=DM_POSTION(i),DM_POSTION(i+1)-1
                    Dynamic_Result_cal.Pool_Volume(i,k)=Dynamic_Result_cal.Pool_Volume(i,k)+ElE_Volume(j,mm)
                END DO
            END DO
            !Dynamic_Result_cal.Pool_Volume(Controlgate_NUM+1,k)=Dynamic_Result_cal.Pool_Volume(Controlgate_NUM+1,k)+ElE_Volume(c_element,mm)

            !***蓄量总和*******
            Dynamic_Result_cal.Total_Volume(k)=0.0
            DO i=1,Com_element_cal.Controlgate_NUM+1
                Dynamic_Result_cal.Total_Volume(k)=Dynamic_Result_cal.Total_Volume(k)+Dynamic_Result_cal.Pool_Volume(i,k)
            END DO

            GateBottom=Com_element_cal.Cogate_cl(1).InvertElev
            GateWidth=Com_element_cal.Cogate_cl(1).ParallelNum*Com_element_cal.Cogate_cl(1).SingleWidth
            BeginBottom=Com_element_cal.Ca_cl(1).InInvertElev
            !算一下湍河节制闸前流速
            v_eachTime(k)=1.2*QQ(firstgate_Index,mm)/(ZZ(firstgate_Index,mm)-GateBottom)/GateWidth!系数用以平衡壅水处流速下降
            !算一下刁河闸后水位与湍河闸前水位的平均值
            h_eachTime(k)=((ZZ(1,mm)-BeginBottom)+(ZZ(firstgate_Index,mm)-GateBottom))/2
            !!47为湍河闸前，严陵河为76
            !v_eachTime(k)=1.2*QQ(76,mm)/(ZZ(76,mm)-138.53)/26!节制闸前流速，27底宽，139.25底高程，当矩形算了，稍微乘个系数平衡壅水处流速较慢
            !!算一下刁河闸后水位与严陵河闸前水位的平均值
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
    !污染带长度
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


    !///////释放空间
    DEALLOCATE(h_new)
    !DEALLOCATE(ZZ)
    !DEALLOCATE(QQ)
    DEALLOCATE(ZZOLD)
    DEALLOCATE(QQOLD)
    DEALLOCATE(SG_POSTION)
    END















    !最后的测试
    SUBROUTINE Dynamic_Cal2(Boundary_temp,pollu_Lenth,arrive_time,Endtime)  !非恒定流计算
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
    INTEGER,ALLOCATABLE::DM_POSTION(:)  !节制闸+首尾断面位置
    INTEGER,ALLOCATABLE::SG_POSTION(:)  !节制闸位置
    REAL::m_coeff    !末端闸门过闸流量系数计算
    REAL::PPQ, VVQ  !水位~流量关系断面系数
    INTEGER::Canal_NUM
    INTEGER::c_element    !计算调整空间步长后参与计算元件的个数

    INTEGER::i,j,k,mm,kk,LLL
    real::c11,c22,c21,hh
    REAL::zzaverage,QQaverage
    INTEGER::si_n   !优化存储时的循环变量
    REAL::random    !随机数
    DOUBLE PRECISION r0
    Integer::pos_sg,DM_sg !记录节制闸位置
    INTEGER::ele_sg    !记录各元件的位置（除了桥墩）
    !稳态各断面水深数组（新）
    REAL,ALLOCATABLE::h_new(:)

    !/////////////非恒定流计算
    !未知量存储
    REAL,ALLOCATABLE::ZZ(:,:)
    REAL,ALLOCATABLE::QQ(:,:)
    REAL,ALLOCATABLE::Pool_Volume(:,:)

    !上一步的一直量存储
    REAL,ALLOCATABLE::ZZOLD(:,:)
    REAL,ALLOCATABLE::QQOLD(:,:)
    REAL,ALLOCATABLE::ELE_Volume(:,:)
    !追赶系数
    REAL,ALLOCATABLE::SS(:)
    REAL,ALLOCATABLE::TT(:)



    Real::GateBottom,GateWidth,BeginBottom

    real::pollu_Lenth
    real::arrive_time
    real::t       !时间变量,单位：min
    real::D  !污染源距分水口500米处距离
    real::S  !计算得到的去染污前端距污染源距离
    real::differ    !误差
    real::D1,D2,w    !迭代计算得出的D
    REAL::v,m,B,h,JJ,g   !系数
    real,allocatable::v_eachTime(:)!每次计算的速度
    real,allocatable::h_eachTime(:)!每次计算的h

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
    m=1.0    !自由倍数取1
    B=26.5
    h=6.25
    JJ=0.00001
    g=pollut_weight !污染物质量。该例子为：10t苯酚水，计算单位:g,后续可手动更换污染物质量
    D=pollut_ConveyDistance    !单位:m



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
    !未知量存储
    ALLOCATE(ZZ(c_element+1,NUMMember))
    ALLOCATE(QQ(c_element+1,NUMMember))
    ALLOCATE(ZZOLD(c_element+1,NUMMember))
    ALLOCATE(QQOLD(c_element+1,NUMMember))
    ALLOCATE(SG_POSTION(Com_element_cal.Controlgate_NUM))
    ALLOCATE(DM_POSTION(Com_element_cal.Controlgate_NUM+2))

    ALLOCATE(Dynamic_Result_cal.Pool_Volume(Com_element_cal.Controlgate_NUM+1,Boundary_2cal.NUM))
    ALLOCATE(ElE_Volume(c_element,NUMMember))
    ALLOCATE(Dynamic_Result_cal.Total_Volume(Boundary_2cal.NUM))
    !///////确定闸门所在的位置//////
    pos_sg=1
    DO j=1,c_element
        IF(Com_element_cal.Ele_Rel(j,2)==6)THEN
            SG_POSTION(pos_sg)=j
            pos_sg=pos_sg+1
        ENDIF
    END DO

    !///////确定闸门+首尾断面所在的位置//////
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
            IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////渠道元件
                ZZ(j,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j,3)).InInvertElev+h_new(j)
                ZZ(j+1,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(j,3)).OutInvertElev+h_new(j+1)
            ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
                ZZ(j,mm)=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j,3)).InInvertElev+h_new(j)
                ZZ(j+1,mm)=Com_element_cal.Insi_cl(Com_element_cal.Ele_Rel(j,3)).OutInvertElev+h_new(j+1)
            ENDIF
        END DO
    END DO
    IF(Boundary_2cal.Upsign==1)THEN
        ZZ(1,1)=Boundary_2cal.UpZ(1)
    ENDIF

    !////////总的计算次数为：Boundary_2cal.NUM,在输入文件中计算过了
    si_n=1
    DO k=1,Boundary_2cal.NUM
        DO mm=1,NUMMember

            DO i=1,c_element,1
                QQOLD(i,mm)=QQ(i,mm)
                if((Com_element_cal.Ele_Rel(i,2)==5).and.( ZZ(i,mm)<Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(i,3)).InInvertElev+0.1)) then
                    ZZOLD(i,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(i,3)).InInvertElev+0.1
                    write(*,*)"==============warning（追赶法计算结果）==============="
                    write(*,*) "第",k,"次计算的时候露底了,将水深赋值为0.1继续计算"
                    write(*,*) "如果连续出现这里的露底问题，会导致结果崩溃。具体原因暂未可知"
                    write(*,*) "可以尝试将计算步长缩短，或者调试一下糙率"
                    write(*,*)"================================================"
                else
                    ZZOLD(i,mm)=ZZ(i,mm)
                end if
            END DO


            !最后一个断面单独处理
            QQOLD(c_element+1,mm)=QQ(c_element+1,mm)
            if((Com_element_cal.Ele_Rel(c_element,2)==5).and.( ZZ(c_element,mm)<Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(c_element,3)).OutInvertElev+0.1)) then
                ZZOLD(c_element+1,mm)=Com_element_cal.Ca_cl(Com_element_cal.Ele_Rel(c_element,3)).OutInvertElev+0.1
            else
                ZZOLD(c_element+1,mm)=ZZ(c_element+1,mm)
            end if







            !#TU-------------------------修改到此结束
            !追赶系数
            ALLOCATE(SS(c_element+1))
            ALLOCATE(TT(c_element+1))
            ALLOCATE(Condition_cal.PP(c_element+1))
            ALLOCATE(Condition_cal.VV(c_element+1))







            isovergate = 0
            !***************生成系数************
            IF(Boundary_2cal.Upsign==0)THEN        !上游边界条件分别选择水位和流量.0表示选择流量边界，1表示上游为水位边界
                !/////////上游边界系数
                Condition_cal.PP(1)=Boundary_2cal.UpQ(k)
                QQ(1,mm)=Boundary_2cal.UpQ(k)
                Condition_cal.VV(1)=0.0
                !/////////中间元件系数
                DO j=1,c_element
                    if((isovergate == 1).and.(k > Endtime))Then
                        QQOLD(j:,m)=150
                        ZZOLD(j,m)=ZZOLD(j,m)+0.6
                        isovergate = 0
                    ENDIF

                    !/////先判断类型再调用对应的函数模块，生成相应的系数。
                    IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////渠道元件
                        Call CanalDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==4)THEN     !////分水口元件
                        Call SideOutFlowDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_2cal.OutDischarge(Com_element_cal.Ele_Rel(j,3),k),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==6)THEN     !////节制闸元件
                        Call ControlGateDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_temp.Opendegree(1,k),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                        isovergate = 1
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==7)THEN  !/////渐变段元件
                        Call AreachangeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
                        Call InvertedsiphonDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==11)THEN !/////桥墩元件
                        Call BridgeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==12)THEN    !/////泵站元件
                        Call PumpingstationDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),C_Bladeangle_cl(Com_element_cal.Ele_Rel(j,3)).Each_PS(:,k),Com_element_cal.PumpST_cl(Com_element_cal.Ele_Rel(j,3)).UNIT_NUM,QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ENDIF
                END DO
                !///////////下游边界条件
                IF(Boundary_2cal.Downsign==1)THEN   !///////////水位边界1
                    ZZ(c_element+1,mm)=Boundary_2cal.DownZ(k)
                ELSEIF(Boundary_2cal.Downsign==2)THEN   !///////////水位~流量关系2    此处的下游边界均使用的DownZ，可能会有问题
                    PPQ=m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))-1.0/(2.0*m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))
                    VVQ=-1.0/(2.0*m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev)))
                    ZZ(c_element+1,mm)=(Condition_cal.PP(c_element+1)-PPQ)/(Condition_cal.VV(c_element+1)-VVQ)
                ENDIF

                !*******回代求解水动力过程********！
                DO j=c_element,1,-1
                    ZZ(j,mm)=SS(j+1)-TT(j+1)*ZZ(j+1,mm)
                    QQ(j+1,mm)=Condition_cal.PP(j+1)-Condition_cal.VV(j+1)*ZZ(j+1,mm)
                END DO
            ELSE    !上游为水位边界Boundary_2cal.Upsign==1
                !/////////上游边界系数
                Condition_cal.PP(1)=Boundary_2cal.UpZ(k)
                ZZ(1,mm)=Boundary_2cal.UpZ(k)
                Condition_cal.VV(1)=0.0
                !/////////中间元件系数
                DO j=1,c_element
                    !/////先判断类型再调用对应的函数模块，生成相应的系数。
                    IF(Com_element_cal.Ele_Rel(j,2)==5)THEN  !/////渠道元件
                        Call CanalDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==4)THEN     !////分水口元件
                        Call SideOutFlowDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_2cal.OutDischarge(Com_element_cal.Ele_Rel(j,3),k),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==6)THEN     !////节制闸元件
                        Call ControlGateDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),Boundary_temp.Opendegree(1,k),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                        isovergate = 1
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==7)THEN  !/////渐变段元件
                        Call AreachangeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
                        Call InvertedsiphonDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==11)THEN !/////桥墩元件
                        Call BridgeDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ELSEIF(Com_element_cal.Ele_Rel(j,2)==12)THEN    !/////泵站元件
                        Call PumpingstationDynamicCoeff(j,Com_element_cal.Ele_Rel(j,3),C_Bladeangle_cl(Com_element_cal.Ele_Rel(j,3)).Each_PS(:,k),Com_element_cal.PumpST_cl(Com_element_cal.Ele_Rel(j,3)).UNIT_NUM,QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j+1),TT(j+1),Condition_cal.PP(j+1),Condition_cal.VV(j+1),ElE_Volume(j,mm),Condition_cal)
                    ENDIF
                END DO
                !///////////下游边界条件
                IF(Boundary_2cal.Downsign==0)THEN   !///////////流量边界0
                    QQ(c_element+1,mm)=Boundary_2cal.DownQ(k)
                ELSEIF(Boundary_2cal.Downsign==1)THEN   !///////////水位边界1
                    QQ(c_element+1,mm)=(Condition_cal.PP(c_element+1)-Boundary_2cal.DownZ(k))/Condition_cal.VV(c_element+1)
                ELSEIF(Boundary_2cal.Downsign==2)THEN   !///////////水位~流量关系边界2
                    PPQ=m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))-1.0/(2.0*m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev))
                    VVQ=-1.0/(2.0*m_coeff*5*Boundary_2cal.DownQ(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Com_element_main.Ca_cl(Canal_NUM).OutInvertElev)))
                    QQ(c_element+1,mm)=(PPQ-VVQ*Condition_cal.PP(c_element+1))/(1-VVQ*Condition_cal.VV(c_element+1))
                ENDIF
                !*******回代求解水动力过程********！
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

            !***各渠池求总和*******
            DO i=1,Com_element_cal.Controlgate_NUM+1
                Dynamic_Result_cal.Pool_Volume(i,k)=0.0
                DO j=DM_POSTION(i),DM_POSTION(i+1)-1
                    Dynamic_Result_cal.Pool_Volume(i,k)=Dynamic_Result_cal.Pool_Volume(i,k)+ElE_Volume(j,mm)
                END DO
            END DO
            !Dynamic_Result_cal.Pool_Volume(Controlgate_NUM+1,k)=Dynamic_Result_cal.Pool_Volume(Controlgate_NUM+1,k)+ElE_Volume(c_element,mm)

            !***蓄量总和*******
            Dynamic_Result_cal.Total_Volume(k)=0.0
            DO i=1,Com_element_cal.Controlgate_NUM+1
                Dynamic_Result_cal.Total_Volume(k)=Dynamic_Result_cal.Total_Volume(k)+Dynamic_Result_cal.Pool_Volume(i,k)
            END DO



            GateBottom=Com_element_cal.Cogate_cl(1).InvertElev
            GateWidth=Com_element_cal.Cogate_cl(1).ParallelNum*Com_element_cal.Cogate_cl(1).SingleWidth
            BeginBottom=Com_element_cal.Ca_cl(1).InInvertElev
            v_eachTime(k)=1.2*QQ(firstgate_Index,mm)/(ZZ(firstgate_Index,mm)-GateBottom)/GateWidth!系数用以平衡壅水处流速下降
            h_eachTime(k)=((ZZ(1,mm)-BeginBottom)+(ZZ(firstgate_Index,mm)-GateBottom))/2




            !*******清空存储空间WRITE(24,*)
            DEALLOCATE(SS)
            DEALLOCATE(TT)
            DEALLOCATE(Condition_cal.PP)
            DEALLOCATE(Condition_cal.VV)
        END DO
    END DO

    !////////////////////////////////////////////////////////////////////////////////////////////////


    v=sum(v_eachTime(:))/Boundary_2cal.NUM
    h=sum(h_eachTime(:))/Boundary_2cal.NUM
    !污染带长度
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
    !///////释放空间
    DEALLOCATE(h_new)
    !DEALLOCATE(ZZ)
    !DEALLOCATE(QQ)
    DEALLOCATE(ZZOLD)
    DEALLOCATE(QQOLD)
    DEALLOCATE(SG_POSTION)


    !******************输出非恒定流计算水动力过程******************
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

    
    !!流量过程
    OPEN(UNIT=28,FILE="OUTPUTFILE/DynamicRESULT_Q.txt")
    WRITE(28,'(*(G0.5,:,"   "))')'时间(分钟)\桩号（km）',(Q_char(i),i=1,firstGateIndex)




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

            !WRITE(28,'(*(G0.5,:,","))')(k-1)*Dynamic_Result_cal.DeltaT/REAL(Out_DeltaT),(Dynamic_Result_cal.ZZ(i,k),Dynamic_Result_cal.QQ(i,k),i=1,Com_element_main.Element_NUM-Com_element_main.Bridge_NUM),Dynamic_Result_cal.Zbd(1,k),Dynamic_Result_cal.Qbd(1,k)!常用
            WRITE(28,'(*(G0.5,:,"   "))')(k-1)*Dynamic_Result_cal.DeltaT/REAL(Out_DeltaT)*Boundary_in.Step/60.0,(Dynamic_Result_cal.QQ(i,k),i=1,firstGateIndex)
            mindex=mindex+1
            gatedischarge(mindex)=Dynamic_Result_cal.QQ(firstGateIndex,k)

            if (temp2>Dynamic_Result_cal.ZZ(10,k)) temp2=Dynamic_Result_cal.ZZ(10,k)
        END IF
    END DO
    CLOSE(28)


    !!水位过程
    OPEN(UNIT=28,FILE="OUTPUTFILE/DynamicRESULT_Z.txt")
    if(Output_sign==1)then
        !WRITE(28,'(*(G0.5,:,"  "))')'时间序列',(Z_char(i),Q_char(i),i=1,Com_element_main.Element_NUM-Com_element_main.Bridge_NUM+1)
        WRITE(28,'(*(G0.5,:,"   "))')'时间(分钟)\桩号（km）',(Q_char(i),i=1,firstGateIndex)
    end if
    DO k=1,Dynamic_Result_cal.NUM,1
        IF(MOD((k-1)*Dynamic_Result_cal.DeltaT,Out_DeltaT)==0)THEN
            WRITE(28,'(*(G0.5,:,"   "))')(k-1)*Dynamic_Result_cal.DeltaT/REAL(Out_DeltaT)*Boundary_in.Step/60.0,(Dynamic_Result_cal.ZZ(i,k),i=1,firstGateIndex)
        END IF
    END DO
    CLOSE(28)


    !输出优化后的节制闸开度
    open(unit=71,file="outputfile/selected-operating.txt")
    write(71,*)"时间(min) 下节制闸开度(mm)  下节制闸Q  退水闸Q(m3/s) 总出流Q   上节制闸Q  "
    do i=1,size(Boundary_2in.Opendegree,2)
        WRITE(71,'(*(G0.5,:,"   "))') (i-1)*Boundary_2in.Step/60,int(Boundary_2in.Opendegree(1,i)*1000),gatedischarge(i),Boundary_2in.OutDischarge(outflow_index,i),gatedischarge(i)+Boundary_2in.OutDischarge(outflow_index,i),Boundary_2in.UPQ(i)
    end do
    close(71)

    !OPEN(UNIT=29,FILE="OUTPUTFILE/out-ini.txt")
    !WRITE(29,'(*(G0.5,:,"   "))')(ZZ(i,1),i=1,53)
    !WRITE(29,'(*(G0.5,:,"   "))')(QQ(i,1),i=1,53)
    !CLOSE(29)

    !write(*,*) '退水闸前最低水位为:', temp2-137.642


    END
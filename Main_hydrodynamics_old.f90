    !Subroutine MAIN_hydrodynamic(Bav,D_Target,T_Target,M_Target)
    Subroutine MAIN_hydrodynamic(Bav,D_Target,M_Target)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    !/////一维明渠水流规则断面水动力模拟及参数率定///////
    INTEGER::i,j,k,gg,fpDynamicResult,fpControlGate,fpQControlGate,fpVolume,&
        fpBestDynamicResult,fpBestControlGate,fpBestQControlGate,fpBestVolume
    REAL::A_n,B_n
    REAL::D_Target  !目标函数1：距目标水位的最大偏差最小
    REAL::M_Target  !目标函数2：距目标水位的平均偏差最小
    !REAL::T_Target   !目标函数3：闸门调控次数最少
    REAL::Bav(1:Controlgate_NUM*Controltime)!Controltime输出时间步长 3secs
    !将水动力计算的全局变量均转变为局部变量
    REAL,ALLOCATABLE::GAV_Opendegree(:) !变量传输过程
    REAL,ALLOCATABLE::Deltadegree(:,:)
    REAL,ALLOCATABLE::Opendegree(:,:)   !闸门开度过程
    REAL,ALLOCATABLE::C_Opendegree(:,:) !计算过程使用闸门开度
    REAL,ALLOCATABLE::ST_STATE(:,:)     !存储各节制闸每小时闸前水位、闸后水位
    REAL,ALLOCATABLE::DS_STATE(:)       !存储最下游边界水位
    REAL,ALLOCATABLE::Pool_Volume(:,:)
    REAL::QControlGate((End_Time-Start_Time)/3600+1,2*Controlgate_NUM) !闸门流量过程
    !记录最后的水力过渡过程
    fpDynamicResult=50
    fpBestDynamicResult=150
    fpControlGate=51
    fpBestControlGate=151
    fpQControlGate=52
    fpBestQControlGate=152
    fpVolume=53
    fpBestVolume=153
    OPEN(UNIT = fpDynamicResult, FILE = ('out/Dynamicresult.txt'))
    OPEN(UNIT = fpControlGate, FILE = ('out/ControlGate.txt'))     !各节制闸闸前闸后水位
    OPEN(UNIT = fpQControlGate, FILE = ('out/QControlGate.txt'))   !各节制闸流量
    OPEN(UNIT = fpVolume, FILE = ('out/Volume.txt'))               !各渠段蓄量
    OPEN(UNIT = fpBestDynamicResult, FILE = ('out/BestDynamicresult.txt'))
    OPEN(UNIT = fpBestControlGate, FILE = ('out/BestControlGate.txt'))     !最佳各节制闸闸前闸后水位
    OPEN(UNIT = fpBestQControlGate, FILE = ('out/BestQControlGate.txt'))   !最佳各节制闸流量
    OPEN(UNIT = fpBestVolume, FILE = ('out/BestVolume.txt'))               !各渠段蓄量
    !*****************************************！
    !*****将调控过程传递到单次水动力计算******！
    !*****************************************！
    ALLOCATE(GAV_Opendegree(Controlgate_NUM*Controltime))
    !    write(*,*)D_Target,M_Target
    DO i=1,Controlgate_NUM*Controltime
        GAV_Opendegree(i)=Bav(i)
    END DO
    ALLOCATE(Deltadegree(Controlgate_NUM,Boundary_NUM))
    ALLOCATE(Opendegree(Controlgate_NUM,Boundary_NUM))
    DeltaDegree(:,:)=0.0  !开始赋值为0
    
    !这一段是什么意思？
    !反正输入文件节制闸只需要第一个开度，其余的会自动补齐？
    DO i=1,Controlgate_NUM
        DO j=1,Controltime    !存在6次调整机会
            DeltaDegree(i,4*j-3)=GAV_Opendegree(Controltime*(i-1)+j)
            DeltaDegree(i,4*j-2)=0.0
            DeltaDegree(i,4*j-1)=0.0
            DeltaDegree(i,4*j-0)=0.0
        END DO
    END DO

    Opendegree(:,1)=Opendegree_INIT(:,1)!输入的初始开度
    DO i=1,Controlgate_NUM
        Opendegree(i,1)=Opendegree(i,1)+DeltaDegree(i,1)
        !行坐标是有几个节制闸，列坐标是有几个计算时段
        DO j=2,Boundary_NUM
            Opendegree(i,j)=Opendegree(i,j-1)+DeltaDegree(i,j)
        ENDDO
    ENDDO
    !////输出观察一下
    ! WRITE(*,*)"当前计算粒子开度调控过程！！！"
    !DO j=1,Boundary_NUM
    !WRITE(*,"(13(F7.3,2X))")(Opendegree(i,j),i=1,Controlgate_NUM)
    !ENDDO

    !//////闸控过程加密，应该放入水动力计算中
    IF(Boundary_Step/=DeltaT)THEN
        ALLOCATE(C_Opendegree(Controlgate_NUM,Calculate_NUM))
        k=Boundary_Step/DeltaT
        DO j=1,Calculate_NUM-2
            C_Opendegree(:,j)=Opendegree(:,j/k+1)+(MOD(j,k)-1)*(Opendegree(:,j/k+2)-Opendegree(:,j/k+1))/k
        END DO
        C_Opendegree(:,Calculate_NUM-1)=Opendegree(:,Boundary_NUM)-(Opendegree(:,Boundary_NUM)-Opendegree(:,Boundary_NUM-1))/k
        C_Opendegree(:,Calculate_NUM)=Opendegree(:,Boundary_NUM)
    ELSEIF(Boundary_Step==DeltaT)THEN
        ALLOCATE(C_Opendegree(Controlgate_NUM,Calculate_NUM))
        k=Boundary_Step/DeltaT
        DO j=1,Calculate_NUM
            C_Opendegree(:,j)=Opendegree(:,j)
        END DO
    END IF
    !/////////////////////////////////////////////////////////////////////////!
    !/////////输出观察一下//////////
    !   WRITE(*,*)"开始非恒定流计算！！！"
    !******************************!
    !*********非恒定流计算*********！
    !******************************!
    
    ALLOCATE(ST_STATE((End_Time-Start_Time)/3600+1,2*Controlgate_NUM))   !存储每小时各节制闸的闸前水位、闸后水位，用于计算目标函数
    ALLOCATE(DS_STATE((End_Time-Start_Time)/3600+1))    !存储下游边界水位，用于计算目标函数
    
    ALLOCATE(Pool_Volume(Controlgate_NUM+1,Calculate_NUM))
    !ALLOCATE(QControlGate((End_Time-Start_Time)/3600+1,2*Controlgate_NUM))
    CALL Dynamic_Cal(C_Opendegree,ST_STATE,DS_STATE,Pool_Volume,QControlGate)!除了开度，其余全是输出值
    !/////非恒定流计算


    !******************************!
    !*********目标函数计算**********！
    !******************************!
    CALL Function_Target(GAV_Opendegree,ST_STATE,DS_STATE,D_Target,M_Target)!,QControlGate



    IF(D_Target<=MIN_Target1)THEN
        write(*,*)D_Target,M_Target
        BestST_STATE=ST_STATE
        BestOpendegree=Opendegree
        BestQControlGate=QControlGate
        BestVolume=Pool_Volume
        write(*,*)"已更新最优值"
    END IF

    !*********输出非恒定流结果********！
    If(endflag==1)then
        !////输出水位变化////
        write(fpDynamicResult,*)"----"
        write(fpDynamicResult,*)D_Target,M_Target
        CALL writeTOfile(fpDynamicResult,ST_STATE,(End_Time-Start_Time)/3600+1,2*Controlgate_NUM)
        !////输出闸门变化////
        write(fpControlGate,*)"----"
        write(fpControlGate,*)D_Target,M_Target
        DO j=1,Boundary_NUM
            WRITE(fpControlGate,"(59(F7.3,2X))")(Opendegree(i,j),i=1,Controlgate_NUM)
        ENDDO
        !////输出闸门流量变化////
        write(fpQControlGate,*)"----"
        write(fpQControlGate,*)D_Target,M_Target
        CALL writeTOfile(fpQControlGate,QControlGate,(End_Time-Start_Time)/3600+1,2*Controlgate_NUM)
        !////输出渠池蓄量变化////
        write(fpVolume,*)"----"
        write(fpVolume,*)D_Target,M_Target
        CALL writeTOfile(fpVolume,Pool_Volume,Controlgate_NUM+1,Calculate_NUM)

        IF(bestoutputflag==1)THEN
            bestoutputflag=0
            !////输出最佳水位变化////
            CALL writeTOfile(fpBestDynamicResult,BestST_STATE,(End_Time-Start_Time)/3600+1,2*Controlgate_NUM)
            !////输出最佳闸门变化////
            DO j=1,Boundary_NUM
                WRITE(fpBestControlGate,"(59(F7.3,2X))")(BestOpendegree(i,j),i=1,Controlgate_NUM)
            ENDDO
            !////输出最佳闸门流量变化////
            CALL writeTOfile(fpBestQControlGate,BestQControlGate,(End_Time-Start_Time)/3600+1,2*Controlgate_NUM)
            !////输出最佳闸门流量变化////
            CALL writeTOfile(fpBestVolume,BestVolume,Controlgate_NUM+1,Calculate_NUM)
        ENDIF
    end if


    DEALLOCATE(C_Opendegree)
    DEALLOCATE(ST_STATE)
    DEALLOCATE(DS_STATE)
    DEALLOCATE(GAV_Opendegree)
    DEALLOCATE(Opendegree)
    DEALLOCATE(Deltadegree)
    DEALLOCATE(Dynamic_Result_Zup)
    DEALLOCATE(Dynamic_Result_Zdown)
    DEALLOCATE(Dynamic_Result_Q)
    DEALLOCATE(Dynamic_Result_Zbd)
    DEALLOCATE(Dynamic_Result_Qbd)
    DEALLOCATE(Dynamic_Result_ZZ)
    DEALLOCATE(Dynamic_Result_QQ)
    END
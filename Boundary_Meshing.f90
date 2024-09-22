!SUBROUTINE Boundary_Meshing(Boundary_in,Boundary_cal)
SUBROUTINE Boundary_Meshing
    USE GLOBAL_VARIABLE
    
    !TYPE(Boundary)::Boundary_in,Boundary_cal
    INTEGER::i,j,k
    
    INTERFACE
        !对不同大小的一维数组进行线性插值
        SUBROUTINE Linear_Interpolation(array1,array2)
            REAL,ALLOCATABLE::array1(:),array2(:)
        END SUBROUTINE
        !对不同大小的二维数组进行线性插值
        SUBROUTINE Linear_Interpolation_2d(array1,array2,direction)
            REAL,ALLOCATABLE::array1(:,:),array2(:,:)
            INTEGER::direction
        END SUBROUTINE
    END INTERFACE
    
    Boundary_cal.Start_Time=Boundary_in.Start_Time
    Boundary_cal.End_Time=Boundary_in.End_Time
    Boundary_cal.Upsign=Boundary_in.Upsign
    Boundary_cal.Downsign=Boundary_in.Downsign
    Boundary_cal.Gate_NUM=Boundary_in.Gate_NUM
    Boundary_cal.PumpST_NUM=Boundary_in.PumpST_NUM
    Boundary_cal.SideOutFlow_NUM=Boundary_in.SideOutFlow_NUM
    !边界条件的插值计算
    CALL Boundary_Auto_Update(Boundary_cal,"NUM")
    CALL Boundary_Auto_Update(Boundary_cal,"Gate_NUM")
    CALL Boundary_Auto_Update(Boundary_cal,"PumpST_NUM")
    CALL Boundary_Auto_Update(Boundary_cal,"SideOutFlow_NUM")
    CALL Linear_Interpolation(Boundary_in.UpQ,Boundary_cal.UpQ)
    CALL Linear_Interpolation(Boundary_in.UpZ,Boundary_cal.UpZ)
    CALL Linear_Interpolation(Boundary_in.DownQ,Boundary_cal.DownQ)
    CALL Linear_Interpolation(Boundary_in.DownZ,Boundary_cal.DownZ)
    CALL Linear_Interpolation_2d(Boundary_in.Opendegree,Boundary_cal.Opendegree,2)
    !CALL Linear_Interpolation_2d(Boundary_in.Bladeangle_cl,Boundary_cal.Bladeangle_cl,2)      !这里还有些问题
    CALL Linear_Interpolation_2d(Boundary_in.OutDischarge,Boundary_cal.OutDischarge,2)


    Calculate_NUM=Boundary_cal.NUM

    
   
        
    !调控过程加密处理
        IF(Boundary_in.Step/=Boundary_cal.Step)THEN
            ALLOCATE(C_Bladeangle_cl(PumpST_NUM))
            DO i=1,PumpST_NUM
                ALLOCATE(C_Bladeangle_cl(i).Each_PS(PumpST_cl(i).UNIT_NUM,Boundary_cal.NUM))
            END DO
            k=Boundary_in.Step/Boundary_cal.Step
            DO j=1,Boundary_cal.NUM-2
                DO i=1,PumpST_NUM
                    C_Bladeangle_cl(i).Each_PS(:,j)=Bladeangle_cl(i).Each_PS(:,j/k+1)
                END DO
            END DO
        ELSEIF(Boundary_in.Step==Boundary_cal.Step)THEN
            DO i=1,PumpST_NUM
                ALLOCATE(C_Bladeangle_cl(i).Each_PS(PumpST_cl(i).UNIT_NUM,Boundary_cal.NUM))
            END DO
            k=Boundary_in.Step/Boundary_cal.Step
            DO j=1,Boundary_cal.NUM
                DO i=1,PumpST_NUM
                    C_Bladeangle_cl(i).Each_PS(:,j)=Bladeangle_cl(i).Each_PS(:,j)
                END DO
            END DO
        END IF
    
    END
    
    
    
    
    
    SUBROUTINE Boundary_Meshing_fun(Boundary_inter,Boundary_calcu)

    USE GLOBAL_VARIABLE
    
    TYPE(Boundary)::Boundary_inter,Boundary_calcu
    INTEGER::i,j,k
    
    INTERFACE
        !对不同大小的一维数组进行线性插值
        SUBROUTINE Linear_Interpolation(array1,array2)
            REAL,ALLOCATABLE::array1(:),array2(:)
        END SUBROUTINE
        !对不同大小的二维数组进行线性插值
        SUBROUTINE Linear_Interpolation_2d(array1,array2,direction)
            REAL,ALLOCATABLE::array1(:,:),array2(:,:)
            INTEGER::direction
        END SUBROUTINE
    END INTERFACE
    
    Boundary_calcu.Start_Time=Boundary_inter.Start_Time
    Boundary_calcu.End_Time=Boundary_inter.End_Time
    Boundary_calcu.Upsign=Boundary_inter.Upsign
    Boundary_calcu.Downsign=Boundary_inter.Downsign
    Boundary_calcu.Gate_NUM=Boundary_inter.Gate_NUM
    Boundary_calcu.PumpST_NUM=Boundary_inter.PumpST_NUM
    Boundary_calcu.SideOutFlow_NUM=Boundary_inter.SideOutFlow_NUM
    !边界条件的插值计算
    CALL Boundary_Auto_Update(Boundary_calcu,"NUM")
    CALL Boundary_Auto_Update(Boundary_calcu,"Gate_NUM")
    CALL Boundary_Auto_Update(Boundary_calcu,"PumpST_NUM")
    CALL Boundary_Auto_Update(Boundary_calcu,"SideOutFlow_NUM")
    CALL Linear_Interpolation(Boundary_inter.UpQ,Boundary_calcu.UpQ)
    CALL Linear_Interpolation(Boundary_inter.UpZ,Boundary_calcu.UpZ)
    CALL Linear_Interpolation(Boundary_inter.DownQ,Boundary_calcu.DownQ)
    CALL Linear_Interpolation(Boundary_inter.DownZ,Boundary_calcu.DownZ)
    CALL Linear_Interpolation_2d(Boundary_inter.Opendegree,Boundary_calcu.Opendegree,2)
    !CALL Linear_Interpolation_2d(Boundary_in.Bladeangle_cl,Boundary_calcu.Bladeangle_cl,2)      !这里还有些问题
    CALL Linear_Interpolation_2d(Boundary_inter.OutDischarge,Boundary_calcu.OutDischarge,2)


    Calculate_NUM=Boundary_calcu.NUM

    
   
        
    !调控过程加密处理
        !IF(Boundary_inter.Step/=Boundary_calcu.Step)THEN
        !    ALLOCATE(C_Bladeangle_cl(PumpST_NUM))
        !    DO i=1,PumpST_NUM
        !        ALLOCATE(C_Bladeangle_cl(i).Each_PS(PumpST_cl(i).UNIT_NUM,Boundary_calcu.NUM))
        !    END DO
        !    k=Boundary_inter.Step/Boundary_calcu.Step
        !    DO j=1,Boundary_calcu.NUM-2
        !        DO i=1,PumpST_NUM
        !            C_Bladeangle_cl(i).Each_PS(:,j)=Bladeangle_cl(i).Each_PS(:,j/k+1)
        !        END DO
        !    END DO
        !ELSEIF(Boundary_inter.Step==Boundary_calcu.Step)THEN
        !    DO i=1,PumpST_NUM
        !        ALLOCATE(C_Bladeangle_cl(i).Each_PS(PumpST_cl(i).UNIT_NUM,Boundary_calcu.NUM))
        !    END DO
        !    k=Boundary_inter.Step/Boundary_calcu.Step
        !    DO j=1,Boundary_calcu.NUM
        !        DO i=1,PumpST_NUM
        !            C_Bladeangle_cl(i).Each_PS(:,j)=Bladeangle_cl(i).Each_PS(:,j)
        !        END DO
        !    END DO
        !END IF
    
END
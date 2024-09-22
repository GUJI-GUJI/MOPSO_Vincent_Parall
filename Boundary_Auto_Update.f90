SUBROUTINE Boundary_Auto_Update(Boundary_old,name)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    
    TYPE(Boundary)::Boundary_old
    CHARACTER(len=*)::name
    
    !更新边界数量并更新数组长度
    IF(TRIM(name)=="NUM")THEN
        Boundary_old.NUM=(Boundary_old.End_Time-Boundary_old.Start_Time)/Boundary_old.Step+1
        IF(SIZE(Boundary_old.UpQ)/=Boundary_old.NUM)THEN
            IF(ALLOCATED(Boundary_old.UpQ))THEN
                DEALLOCATE(Boundary_old.UpQ)
                DEALLOCATE(Boundary_old.UpZ)
                DEALLOCATE(Boundary_old.DownQ)
                DEALLOCATE(Boundary_old.DownZ)
            END IF
            ALLOCATE(Boundary_old.UpQ(Boundary_old.NUM))
            ALLOCATE(Boundary_old.UpZ(Boundary_old.NUM))
            ALLOCATE(Boundary_old.DownQ(Boundary_old.NUM))
            ALLOCATE(Boundary_old.DownZ(Boundary_old.NUM))
        END IF
    ELSE IF(TRIM(name)=="Gate_NUM")THEN
        IF(SIZE(Boundary_old.Opendegree(:,1))/=Boundary_old.Gate_NUM.or.&
            SIZE(Boundary_old.Opendegree(1,:))/=Boundary_old.NUM)THEN
            IF(ALLOCATED(Boundary_old.Opendegree))THEN
                DEALLOCATE(Boundary_old.Opendegree)
            END IF
            ALLOCATE(Boundary_old.Opendegree(Boundary_old.Gate_NUM,Boundary_old.NUM))
        END IF
    ELSE IF(TRIM(name)=="SideOutFlow_NUM")THEN
        IF(SIZE(Boundary_old.OutDischarge(:,1))/=Boundary_old.SideOutFlow_NUM.or.&
            SIZE(Boundary_old.OutDischarge(1,:))/=Boundary_old.NUM)THEN
            IF(ALLOCATED(Boundary_old.OutDischarge))THEN
                DEALLOCATE(Boundary_old.OutDischarge)
            END IF
            ALLOCATE(Boundary_old.OutDischarge(Boundary_old.SideOutFlow_NUM,Boundary_old.NUM))
        END IF
    ELSE IF(TRIM(name)=="PumpST_NUM")THEN
        IF(SIZE(Boundary_old.Bladeangle_cl(:))/=Boundary_old.PumpST_NUM)THEN
            IF(ALLOCATED(Boundary_old.Bladeangle_cl))THEN
                DEALLOCATE(Boundary_old.Bladeangle_cl)
            END IF
            ALLOCATE(Boundary_old.Bladeangle_cl(Boundary_old.PumpST_NUM))
        END IF       
    END IF
    
END
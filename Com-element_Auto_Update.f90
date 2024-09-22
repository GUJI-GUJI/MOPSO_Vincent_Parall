SUBROUTINE Com_element_Auto_Update(Com_element_old,name)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    
    TYPE(Com_element)::Com_element_old
    CHARACTER(len=*)::name
    
    IF(TRIM(name)=="NUM")THEN
        IF(SIZE(Com_element_old.Ele_Rel(:,1))/=Com_element_old.Element_NUM)THEN
            IF(ALLOCATED(Com_element_old.Ele_Rel))THEN
                DEALLOCATE(Com_element_old.Ele_Rel)
            END IF
            ALLOCATE(Com_element_old.Ele_Rel(Com_element_old.Element_NUM,3))
        END IF
        IF(SIZE(Com_element_old.Ca_cl)/=Com_element_old.Canal_NUM)THEN
            IF(ALLOCATED(Com_element_old.Ca_cl))THEN
                DEALLOCATE(Com_element_old.Ca_cl)
            END IF
            ALLOCATE(Com_element_old.Ca_cl(Com_element_old.Canal_NUM))
        END IF
        IF(SIZE(Com_element_old.Arch_cl)/=Com_element_old.Areachange_NUM)THEN
            IF(ALLOCATED(Com_element_old.Arch_cl))THEN
                DEALLOCATE(Com_element_old.Arch_cl)
            END IF
            ALLOCATE(Com_element_old.Arch_cl(Com_element_old.Areachange_NUM))
        END IF
        IF(SIZE(Com_element_old.Insi_cl)/=Com_element_old.Invertedsiphon_NUM)THEN
            IF(ALLOCATED(Com_element_old.Insi_cl))THEN
                DEALLOCATE(Com_element_old.Insi_cl)
            END IF
            ALLOCATE(Com_element_old.Insi_cl(Com_element_old.Invertedsiphon_NUM))
        END IF
        IF(SIZE(Com_element_old.Cogate_cl)/=Com_element_old.Controlgate_NUM)THEN
            IF(ALLOCATED(Com_element_old.Cogate_cl))THEN
                DEALLOCATE(Com_element_old.Cogate_cl)
            END IF
            ALLOCATE(Com_element_old.Cogate_cl(Com_element_old.Controlgate_NUM))
        END IF
        IF(SIZE(Com_element_old.SideOF_cl)/=Com_element_old.SideOutFlow_NUM)THEN
            IF(ALLOCATED(Com_element_old.SideOF_cl))THEN
                DEALLOCATE(Com_element_old.SideOF_cl)
            END IF
            ALLOCATE(Com_element_old.SideOF_cl(Com_element_old.SideOutFlow_NUM))
        END IF
        IF(SIZE(Com_element_old.Brid_cl)/=Com_element_old.Bridge_NUM)THEN
            IF(ALLOCATED(Com_element_old.Brid_cl))THEN
                DEALLOCATE(Com_element_old.Brid_cl)
            END IF
            ALLOCATE(Com_element_old.Brid_cl(Com_element_old.Bridge_NUM))
        END IF
        IF(SIZE(Com_element_old.PumpST_cl)/=Com_element_old.PumpST_NUM)THEN
            IF(ALLOCATED(Com_element_old.PumpST_cl))THEN
                DEALLOCATE(Com_element_old.PumpST_cl)
            END IF
            ALLOCATE(Com_element_old.PumpST_cl(Com_element_old.PumpST_NUM))
        END IF
        IF(SIZE(Com_element_old.n(:))/=Com_element_old.Element_NUM)THEN
            IF(ALLOCATED(Com_element_old.n))THEN
                DEALLOCATE(Com_element_old.n)
            END IF
            ALLOCATE(Com_element_old.n(Com_element_old.Element_NUM))
        END IF
        
    ELSE IF(TRIM(name)=="Pump_type_NUM")THEN
        IF(SIZE(Com_element_old.Pump_cl)/=Com_element_old.Pump_type_NUM)THEN
            IF(ALLOCATED(Com_element_old.Pump_cl))THEN
                DEALLOCATE(Com_element_old.Pump_cl)
            END IF
            ALLOCATE(Com_element_old.Pump_cl(Com_element_old.Pump_type_NUM))
        END IF
    END IF
    
    
END
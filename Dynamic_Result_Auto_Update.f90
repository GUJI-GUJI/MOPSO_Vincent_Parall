SUBROUTINE Dynamic_Result_Auto_Update(Dynamic_Result_old,Com_element_obj,Boundary_obj,name)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    
    TYPE(Dynamic_Result)::Dynamic_Result_old
    TYPE(Com_element)::Com_element_obj
    TYPE(Boundary)::Boundary_obj
    CHARACTER(len=*)::name

    IF(TRIM(name)=="NUM")THEN
        IF(SIZE(Dynamic_Result_old.Zup(:,1))/=Com_element_obj.Controlgate_NUM.or.&
            SIZE(Dynamic_Result_old.Zup(1,:))/=Boundary_obj.NUM)THEN
            IF(ALLOCATED(Dynamic_Result_old.Zup))THEN
                DEALLOCATE(Dynamic_Result_old.Zup)
                DEALLOCATE(Dynamic_Result_old.Zdown)
                DEALLOCATE(Dynamic_Result_old.Q)
                DEALLOCATE(Dynamic_Result_old.Zbd)
                DEALLOCATE(Dynamic_Result_old.Qbd)
            END IF
            ALLOCATE(Dynamic_Result_old.Zup(Com_element_obj.Controlgate_NUM,Boundary_obj.NUM))
            ALLOCATE(Dynamic_Result_old.Zdown(Com_element_obj.Controlgate_NUM,Boundary_obj.NUM))
            ALLOCATE(Dynamic_Result_old.Q(Com_element_obj.Controlgate_NUM,Boundary_obj.NUM))
            ALLOCATE(Dynamic_Result_old.Zbd(1,Boundary_obj.NUM))
            ALLOCATE(Dynamic_Result_old.Qbd(1,Boundary_obj.NUM))
        END IF
        
        IF(SIZE(Dynamic_Result_old.ZZ(:,1))/=Com_element_obj.Element_NUM-Com_element_obj.Bridge_NUM.or.&
            SIZE(Dynamic_Result_old.ZZ(1,:))/=Boundary_obj.NUM)THEN
            IF(ALLOCATED(Dynamic_Result_old.ZZ))THEN
                DEALLOCATE(Dynamic_Result_old.ZZ)
                DEALLOCATE(Dynamic_Result_old.QQ)
            END IF
            ALLOCATE(Dynamic_Result_old.ZZ(Com_element_obj.Element_NUM-Com_element_obj.Bridge_NUM,Boundary_obj.NUM))
            ALLOCATE(Dynamic_Result_old.QQ(Com_element_obj.Element_NUM-Com_element_obj.Bridge_NUM,Boundary_obj.NUM))
        END IF
    END IF
    
END
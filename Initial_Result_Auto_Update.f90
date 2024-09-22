SUBROUTINE Initial_Result_Auto_Update(Initial_Result_old,Com_element_obj,name)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    
    TYPE(Initial_Result)::Initial_Result_old
    TYPE(Com_element)::Com_element_obj
    CHARACTER(len=*)::name
    
    IF(TRIM(name)=="NUM")THEN
        IF(SIZE(Initial_Result_old.h)/=Com_element_obj.Element_NUM+1)THEN
            IF(ALLOCATED(Initial_Result_old.h))THEN
                DEALLOCATE(Initial_Result_old.h)
                DEALLOCATE(Initial_Result_old.Q)
                DEALLOCATE(Initial_Result_old.Z)
                DEALLOCATE(Initial_Result_old.vc)
                DEALLOCATE(Initial_Result_old.L)
                DEALLOCATE(Initial_Result_old.vol)
            END IF
            ALLOCATE(Initial_Result_old.h(Com_element_obj.Element_NUM+1))
            ALLOCATE(Initial_Result_old.Q(Com_element_obj.Element_NUM+1))
            ALLOCATE(Initial_Result_old.Z(Com_element_obj.Element_NUM+1))
            ALLOCATE(Initial_Result_old.vc(Com_element_obj.Element_NUM+1))
            ALLOCATE(Initial_Result_old.L(Com_element_obj.Element_NUM+1))
            ALLOCATE(Initial_Result_old.vol(Com_element_obj.Element_NUM+1))
        END IF
        IF(SIZE(Initial_Result_old.Pool_Volume)/=Com_element_obj.Controlgate_NUM+1)THEN
            IF(ALLOCATED(Initial_Result_old.Pool_Volume))THEN
                DEALLOCATE(Initial_Result_old.Pool_Volume)
            END IF
            ALLOCATE(Initial_Result_old.Pool_Volume(Com_element_obj.Controlgate_NUM+1))
        END IF
    END IF
    
END
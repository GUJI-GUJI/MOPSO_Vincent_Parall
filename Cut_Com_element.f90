SUBROUTINE Cut_Com_element(Com_element_in,Com_element_out,Start_Element,End_Element)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!用于从长距离的模型中截取其中一部分!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Com_element_in:长距离的模型
    !Com_element_out:截取得到的单闸模型
    !Start_Element:首端元件（包含）
    !End_Element:末端元件（包含）
    
    !使用条件：
    !1.元件关系使用Ele_Rel描述；
    !2.仅包含渠道、渐变段、倒虹吸、桥墩、闸门、分水口和泵站元件
    
    !适用HOMS版本：Beta2.0.1~Beta2.0.1
    
    !23/2/17 李汉元
    
    TYPE(Com_element)::Com_element_in,Com_element_out
    INTEGER::i,Start_Element,End_Element
    INTEGER::Ca_cl_flag,Arch_cl_flag,Insi_cl_flag,Brid_cl_flag,Cogate_cl_flag,SideOF_cl_flag,PumpST_cl_flag
    
    Ca_cl_flag=0
    Arch_cl_flag=0
    Insi_cl_flag=0
    Brid_cl_flag=0
    Cogate_cl_flag=0
    SideOF_cl_flag=0
    PumpST_cl_flag=0
    
    !统计各元件数量
    Com_element_out.Element_NUM=END_Element-Start_Element+1
    Com_element_out.Canal_NUM=0
    Com_element_out.Areachange_NUM=0
    Com_element_out.Invertedsiphon_NUM=0
    Com_element_out.Bridge_NUM=0
    Com_element_out.Controlgate_NUM=0
    Com_element_out.SideOutFlow_NUM=0
    Com_element_out.PumpST_NUM=0
    DO i=Start_Element,END_Element
        IF(Com_element_in.Ele_Rel(i,2)==5)  Com_element_out.Canal_NUM=Com_element_out.Canal_NUM+1
        IF(Com_element_in.Ele_Rel(i,2)==7)  Com_element_out.Areachange_NUM=Com_element_out.Areachange_NUM+1
        IF(Com_element_in.Ele_Rel(i,2)==11) Com_element_out.Bridge_NUM=Com_element_out.Bridge_NUM+1
        IF(Com_element_in.Ele_Rel(i,2)==10) Com_element_out.Invertedsiphon_NUM=Com_element_out.Invertedsiphon_NUM+1
        IF(Com_element_in.Ele_Rel(i,2)==6)  Com_element_out.Controlgate_NUM=Com_element_out.Controlgate_NUM+1
        IF(Com_element_in.Ele_Rel(i,2)==4)  Com_element_out.SideOutFlow_NUM=Com_element_out.SideOutFlow_NUM+1
        IF(Com_element_in.Ele_Rel(i,2)==12) Com_element_out.PumpST_NUM=Com_element_out.PumpST_NUM+1
    ENDDO
    CALL Com_element_Auto_Update(Com_element_out,"NUM")
    
    !将原有数据放进新的com_element里
    Com_element_out.Ele_Rel=0
    DO i=Start_Element,End_Element
        Com_element_out.Ele_Rel(i-Start_Element+1,1)=i-Start_Element+1
        Com_element_out.Ele_Rel(i-Start_Element+1,2)=Com_element_in.Ele_Rel(i,2)
        IF(Com_element_in.Ele_Rel(i,2)==5)THEN
            Ca_cl_flag=Ca_cl_flag+1
            Com_element_out.Ele_Rel(i-Start_Element+1,3)=Ca_cl_flag
            Com_element_out.Ca_cl(Ca_cl_flag)=Com_element_in.Ca_cl(Com_element_in.Ele_Rel(i,3))
        ENDIF
        IF(Com_element_in.Ele_Rel(i,2)==7)THEN
            Arch_cl_flag=Arch_cl_flag+1
            Com_element_out.Ele_Rel(i-Start_Element+1,3)=Arch_cl_flag
            Com_element_out.Arch_cl(Arch_cl_flag)=Com_element_in.Arch_cl(Com_element_in.Ele_Rel(i,3))
        ENDIF
        IF(Com_element_in.Ele_Rel(i,2)==11)THEN
            Brid_cl_flag=Brid_cl_flag+1
            Com_element_out.Ele_Rel(i-Start_Element+1,3)=Brid_cl_flag
            Com_element_out.Brid_cl(Brid_cl_flag)=Com_element_in.Brid_cl(Com_element_in.Ele_Rel(i,3))
        ENDIF
        IF(Com_element_in.Ele_Rel(i,2)==10)THEN
            Insi_cl_flag=Insi_cl_flag+1
            Com_element_out.Ele_Rel(i-Start_Element+1,3)=Insi_cl_flag
            Com_element_out.Insi_cl(Insi_cl_flag)=Com_element_in.Insi_cl(Com_element_in.Ele_Rel(i,3))
        ENDIF
        IF(Com_element_in.Ele_Rel(i,2)==6)THEN
            Cogate_cl_flag=Cogate_cl_flag+1
            Com_element_out.Ele_Rel(i-Start_Element+1,3)=Cogate_cl_flag
            Com_element_out.Cogate_cl(Cogate_cl_flag)=Com_element_in.Cogate_cl(Com_element_in.Ele_Rel(i,3))
        ENDIF
        IF(Com_element_in.Ele_Rel(i,2)==4)THEN
            SideOF_cl_flag=SideOF_cl_flag+1
            Com_element_out.Ele_Rel(i-Start_Element+1,3)=SideOF_cl_flag
            Com_element_out.SideOF_cl(SideOF_cl_flag)=Com_element_in.SideOF_cl(Com_element_in.Ele_Rel(i,3))
        ENDIF
        IF(Com_element_in.Ele_Rel(i,2)==12)THEN
            PumpST_cl_flag=PumpST_cl_flag+1
            Com_element_out.Ele_Rel(i-Start_Element+1,3)=PumpST_cl_flag
            Com_element_out.PumpST_cl(PumpST_cl_flag)=Com_element_in.PumpST_cl(Com_element_in.Ele_Rel(i,3))
        ENDIF
    ENDDO
    
END SUBROUTINE
    
SUBROUTINE Cut_Com_element_Single_Controlgate(Com_element_in,Com_element_out,Controlgate_Serial)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!用于从长距离的模型中截取单闸模型!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Com_element_in:长距离的模型
    !Controlgate_Serial:需要截取的闸门编号
    !Com_element_out:截取得到的单闸模型
    
    !使用条件：
    !同Cut_Com_element
    
    !适用HOMS版本：Beta2.0.1~Beta2.0.1
    
    !23/2/17 李汉元
    
    TYPE(Com_element)::Com_element_in,Com_element_out
    INTEGER::i,Start_Element,End_Element,Controlgate_Serial,found_flag
    
    Start_Element=1
    End_Element=Com_element_in.Element_NUM
    found_flag=0
    DO i=1,Com_element_in.Element_NUM
        IF(Com_element_in.Ele_Rel(i,2)==6.and.Com_element_in.Ele_Rel(i,3)==Controlgate_Serial-1)THEN    !寻找目标闸门的上一个闸门
            Start_Element=i+1
        END IF
        IF(Com_element_in.Ele_Rel(i,2)==6.and.Com_element_in.Ele_Rel(i,3)==Controlgate_Serial)THEN      !寻找目标闸门
            found_flag=1
        END IF
        IF(Com_element_in.Ele_Rel(i,2)==6.and.Com_element_in.Ele_Rel(i,3)==Controlgate_Serial+1)THEN    !寻找目标闸门的下一个闸门
            End_Element=i-1
        END IF
    END DO
    IF(found_flag==0)THEN
        WRITE(*,*)"Target control gate not found when cutting Com_element",Com_element_in.name          !未找到目标闸门的警告
        pause
    END IF
    
    CALL Cut_Com_element(Com_element_in,Com_element_out,Start_Element,End_Element)
    
END SUBROUTINE
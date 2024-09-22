SUBROUTINE Model_Meshing
    USE GLOBAL_VARIABLE
    
    !TYPE(Com_element)::Com_element_main,Com_element_cal
    !TYPE(Initial_Result)::Initial_Result_main
    INTEGER::i,j,k
    
    INTEGER::Canal_New_NUM  !ʵ�ʲ���������εĸ���
    INTEGER::c_element    !��������ռ䲽����������Ԫ���ĸ���
    
    !******************************!
    !******�Ǻ㶨������Ԥ����******��
    !******************************!
    !����Ǻ㶨������ռ估����Ǻ㶨��������µ�Ԫ����Ϣ
    !//////����ϵ������ռ䣬֮�������ֵ���������������ģ��ĵ���
    !//////���㷽�̸�����ȡ����n(i)
    Element_NUM=Com_element_main.Element_NUM
    c_element=0
    DO i=1,Element_NUM,1
        c_element=c_element+Com_element_main.n(i)
    END DO
    
    !//////////���µ�Ԫ�����б�ţ�����ȷ������Ԫ������ػ���������������������Ӱ�쵽�����
    !��Ҫ�Ƕ��������½��оֲ���ţ�����Ԫ�������������±��
    ALLOCATE(Com_element_cal.Ele_Rel(c_element,3))
    i=Element_NUM-Com_element_main.Bridge_NUM
    ALLOCATE(Com_element_cal.ELE_POSTION(i))
    k=0
    gg=1
    DO i=1,Element_NUM
        Com_element_cal.ELE_POSTION(i)=k+1
        DO j=1,Com_element_main.n(i)
            k=k+1
            Com_element_cal.Ele_Rel(k,1)=k
            Com_element_cal.Ele_Rel(k,2)=Com_element_main.Ele_Rel(i,2)
            IF(Com_element_main.Ele_Rel(i,2)==5)THEN
                Com_element_cal.Ele_Rel(k,3)=j-1+gg
            ELSE
                Com_element_cal.Ele_Rel(k,3)=Com_element_main.Ele_Rel(i,3)
            END IF
        END DO
        IF(Com_element_main.Ele_Rel(i,2)==5)THEN
            gg=gg+Com_element_main.n(i)
        END IF
    END DO

    !ȷ���µļ���Ԫ���Ļ���������ʵ����ֻ��ȷ���µ���������Ϣ������Ԫ������ֱ�ӵȼ�
    Com_element_cal.Canal_NUM=gg-1
    ALLOCATE(Com_element_main.n_canal(Com_element_main.Canal_NUM))
    j=0
    DO i=1,Element_NUM
        IF(Com_element_main.Ele_Rel(i,2)==5)THEN
            j=j+1
            Com_element_main.n_canal(j)=Com_element_main.n(i)
        END IF
        !if(Com_element_main.Ele_Rel(i,2)==6)THEN
        !    firstgate_Index = 2*i-1
        !END iF
    END DO
    ALLOCATE(Com_element_cal.Ca_cl(Com_element_cal.Canal_NUM))
    k=0
    DO i=1,Com_element_main.Canal_NUM
        DO j=1,Com_element_main.n_canal(i)
            k=k+1
            Com_element_cal.Ca_cl(k).Serial_NUM=k
            Com_element_cal.Ca_cl(k).Name=Com_element_main.Ca_cl(i).Name
            Com_element_cal.Ca_cl(k).BottomWidth=Com_element_main.Ca_cl(i).BottomWidth
            Com_element_cal.Ca_cl(k).SideSlopes=Com_element_main.Ca_cl(i).SideSlopes
            Com_element_cal.Ca_cl(k).Incoord=Com_element_main.Ca_cl(i).Incoord+(j-1)*(Com_element_main.Ca_cl(i).Outcoord-Com_element_main.Ca_cl(i).Incoord)/Com_element_main.n_canal(i)
            Com_element_cal.Ca_cl(k).Outcoord=Com_element_main.Ca_cl(i).Incoord+j*(Com_element_main.Ca_cl(i).Outcoord-Com_element_main.Ca_cl(i).Incoord)/Com_element_main.n_canal(i)
            Com_element_cal.Ca_cl(k).InInvertElev=Com_element_main.Ca_cl(i).InInvertElev+(j-1)*(Com_element_main.Ca_cl(i).OutInvertElev-Com_element_main.Ca_cl(i).InInvertElev)/Com_element_main.n_canal(i)
            Com_element_cal.Ca_cl(k).OutInvertElev=Com_element_main.Ca_cl(i).InInvertElev+j*(Com_element_main.Ca_cl(i).OutInvertElev-Com_element_main.Ca_cl(i).InInvertElev)/Com_element_main.n_canal(i)
            Com_element_cal.Ca_cl(k).InTopElev=Com_element_main.Ca_cl(i).InTopElev+(j-1)*(Com_element_main.Ca_cl(i).OutTopElev-Com_element_main.Ca_cl(i).InTopElev)/Com_element_main.n_canal(i)
            Com_element_cal.Ca_cl(k).OutTopElev=Com_element_main.Ca_cl(i).InTopElev+j*(Com_element_main.Ca_cl(i).OutTopElev-Com_element_main.Ca_cl(i).InTopElev)/Com_element_main.n_canal(i)
            Com_element_cal.Ca_cl(k).Roughness=Com_element_main.Ca_cl(i).Roughness
        END DO
    END DO
    Com_element_cal.Element_NUM=c_element
    Com_element_cal.Areachange_NUM=Com_element_main.Areachange_NUM
    Com_element_cal.Arch_cl=Com_element_main.Arch_cl
    Com_element_cal.Bridge_NUM=Com_element_main.Bridge_NUM
    Com_element_cal.Brid_cl=Com_element_main.Brid_cl
    Com_element_cal.Controlgate_NUM=Com_element_main.Controlgate_NUM
    Com_element_cal.Cogate_cl=Com_element_main.Cogate_cl
    Com_element_cal.Invertedsiphon_NUM=Com_element_main.Invertedsiphon_NUM
    Com_element_cal.Insi_cl=Com_element_main.Insi_cl
    Com_element_cal.SideOutFlow_NUM=Com_element_main.SideOutFlow_NUM
    Com_element_cal.SideOF_cl=Com_element_main.SideOF_cl
    Com_element_cal.PumpST_NUM=Com_element_main.PumpST_NUM
    Com_element_cal.PumpST_cl=Com_element_main.PumpST_cl
    Com_element_cal.Pump_type_NUM=Com_element_main.Pump_type_NUM
    Com_element_cal.Pump_cl=Com_element_main.Pump_cl
    !*************************************************************************************************************!
    
END
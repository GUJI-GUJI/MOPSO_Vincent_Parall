SUBROUTINE Model_Meshing
    USE GLOBAL_VARIABLE
    
    !TYPE(Com_element)::Com_element_main,Com_element_cal
    !TYPE(Initial_Result)::Initial_Result_main
    INTEGER::i,j,k
    
    INTEGER::Canal_New_NUM  !实际参与计算渠段的个数
    INTEGER::c_element    !计算调整空间步长后参与计算元件的个数
    
    !******************************!
    !******非恒定流计算预处理******！
    !******************************!
    !分配非恒定流计算空间及参与非恒定流计算的新的元件信息
    !//////分配系数矩阵空间，之后逐个赋值，里面包含各个子模块的调用
    !//////计算方程个数，取决于n(i)
    Element_NUM=Com_element_main.Element_NUM
    c_element=0
    DO i=1,Element_NUM,1
        c_element=c_element+Com_element_main.n(i)
    END DO
    
    !//////////对新的元件进行编号，另外确定计算元件的相关基本参数。（计算结果不能影响到输出）
    !主要是对渠道重新进行局部编号，所有元件进行整体重新编号
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

    !确定新的计算元件的基本参数，实际上只用确定新的渠道的信息，其他元件可以直接等价
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
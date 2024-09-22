SUBROUTINE Initial_Cal  !wjb0413
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    
    INTEGER::i,j,k,gg,mm
    !/////初始计算
    !给各断面水位和流量分配存储空间
    ALLOCATE(h(Element_NUM+1))
    ALLOCATE(Q(Element_NUM+1))
    ALLOCATE(Z(Element_NUM+1))
    ALLOCATE(vc(Element_NUM+1))
    ALLOCATE(L(Element_NUM+1))
    ALLOCATE(Vol(Element_NUM+1))
    !给各元件存储个数分配存储空间
    ALLOCATE(nnn(Element_NUM))
	Waterdepth_flag=0
    !/////初始各断面流量计算
    IF(Boundary_sign_UP==0)THEN    !根据上游边界的不同选择流量赋值
        Q(1)=Upstream(1)
    ELSE
        Q(1)=Initial_Discharge
    END IF
    j=0
    DO i=2,Element_NUM
        IF(Ele_Rel(i,2)==4)THEN
            j=j+1
            Q(i)=Q(i-1)-OutDischarge(j,1)
        Else
            Q(i)=Q(i-1)
        END IF
    END DO
    Q(Element_NUM+1)=Q(Element_NUM)

    !//////初始水位的计算，以及各渠段的分段！目的是确定各断面的水位
    !//////水位的计算是从下游往上游推求
    !//////下游水位边界直接赋值水位
    h(Element_NUM+1)=D_elevation
    Waterdepth_flag=Waterdepth_flag+1
    Waterdepth(Waterdepth_flag,1)=h(Element_NUM+1)
    Waterdepth(Waterdepth_flag,2)=Q(Element_NUM+1)
    
        !///////中间各断面水位赋值
        DO j=Element_NUM,1,-1
           !/////先判断类型再调用对应的函数模块，计算各元件入口水深，同时计算各元件的小断面个数。
           IF(Ele_Rel(j,2)==5)THEN  !/////渠道元件
               Call CanalInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ELSEIF(Ele_Rel(j,2)==4)THEN  !/////分数口元件
               Call SideOutFlowInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ELSEIF(Ele_Rel(j,2)==6)THEN  !/////节制闸元件
               Call ControlGateInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ELSEIF(Ele_Rel(j,2)==7)THEN  !/////渐变段元件
               Call AreachangeInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ELSEIF(Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
               Call InvertedsiphonInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
	ELSEIF(Ele_Rel(j,2)==11)THEN !/////桥梁元件
               Call BridgeInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ELSEIF(Ele_Rel(j,2)==12)THEN !/////泵站元件
               Call PumpingstationInitialStageCalc(j,Ele_Rel(j,3),h(j+1),Q(j+1),Q(j),h(j),nnn(j),Vol(j))
           ENDIF 
        END DO
 
    !***************************************!
    !***********输出初始计算结果************!
    !***************************************!
    
    !//////各断面的里程///////////
    DO j=1,Element_NUM,1
        IF(Ele_Rel(j,2)==5)THEN  !/////渠道元件
           L(j)=Ca_cl(Ele_Rel(j,3))%Incoord
           L(j+1)=Ca_cl(Ele_Rel(j,3))%Outcoord
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
           L(j)=Insi_cl(Ele_Rel(j,3))%Incoord
           L(j+1)=Insi_cl(Ele_Rel(j,3))%Outcoord
        ENDIF     
    END DO
    
    !//////各断面的流速///////////
    DO j=1,Element_NUM,1
        IF(Ele_Rel(j,2)==5)THEN  !/////渠道元件
           vc(j)=Q(j)/((Ca_cl(Ele_Rel(j,3))%BottomWidth+Ca_cl(Ele_Rel(j,3))%SideSlopes*h(j))*h(j))
           vc(j+1)=Q(j+1)/((Ca_cl(Ele_Rel(j,3))%BottomWidth+Ca_cl(Ele_Rel(j,3))%SideSlopes*h(j+1))*h(j+1))
        ELSEIF(Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
           vc(j)=Q(j)/(Insi_cl(Ele_Rel(j,3))%ParallelNum*Insi_cl(Ele_Rel(j,3))%BottomWidth*h(j))
           vc(j+1)=Q(j+1)/(Insi_cl(Ele_Rel(j,3))%ParallelNum*Insi_cl(Ele_Rel(j,3))%BottomWidth*h(j+1))
        ENDIF     
    END DO
    
    !//////计算各断面水位，即计算得到的水深加上各断面的底高程
    DO j=1,Element_NUM,1
       IF(Ele_Rel(j,2)==5)THEN  !/////渠道元件
           Z(j)=Ca_cl(Ele_Rel(j,3))%InInvertElev+h(j)
           Z(j+1)=Ca_cl(Ele_Rel(j,3))%OutInvertElev+h(j+1)
       ELSEIF(Ele_Rel(j,2)==10)THEN !/////倒虹吸元件
           Z(j)=Insi_cl(Ele_Rel(j,3))%InInvertElev+h(j)
           Z(j+1)=Insi_cl(Ele_Rel(j,3))%OutInvertElev+h(j+1)
       ENDIF 
    END DO
	Vol(Element_NUM+1)=0
    
    !/////输出初始计算水动力过程//////////
    OPEN(UNIT=20,FILE="OUTPUTFILE/InitialRESULT.txt")
        WRITE(20,"(A10,3X,A10,3X,A7,3X,A7,3X,A7,3X,A7,3X,A7,3X)")"断面序号","里程","水深","流量","流速","水位","总水头"
        DO i=1,Element_NUM+1
            WRITE(20,"(I6,7X,F10.3,3X,F7.3,3X,F7.3,3X,F7.3,3X,F7.3,3X,F7.3,3X)")i,L(i),h(i),Q(i),vc(i),Z(i),Z(i)+vc(i)**2/2.0/GRAV
        END DO
    CLOSE(20)

    !******************************!
    !******非恒定流计算预处理******！
    !******************************!
    !分配非恒定流计算空间及参与非恒定流计算的新的元件信息
    !//////分配系数矩阵空间，之后逐个赋值，里面包含各个子模块的调用
    !//////计算方程个数，取决于n(i)
    c_element=0
    DO i=1,Element_NUM,1
        c_element=c_element+nnn(i)
    END DO
    
    !//////////对新的元件进行编号，另外确定计算元件的相关基本参数。（计算结果不能影响到输出）
    !主要是对渠道重新进行局部编号，所有元件进行整体重新编号
    ALLOCATE(Ele_Rel_New(c_element,3))
    k=0
    gg=1
    DO i=1,Element_NUM
        DO j=1,nnn(i)
            k=k+1
            Ele_Rel_New(k,1)=k
            Ele_Rel_New(k,2)=Ele_Rel(i,2)
            IF(Ele_Rel(i,2)==5)THEN
                Ele_Rel_New(k,3)=j-1+gg
            ELSE
                Ele_Rel_New(k,3)=Ele_Rel(i,3)
            END IF
        END DO
        IF(Ele_Rel(i,2)==5)THEN
            gg=gg+nnn(i)
        END IF
    END DO
!    OPEN(UNIT=21,FILE="TEST.txt")
!        DO i=1,c_element
!            WRITE(21,*)(Ele_Rel_New(i,j),j=1,3)   
!        END DO
!    CLOSE(21)

    !确定新的计算元件的基本参数，实际上只用确定新的渠道的信息，其他元件可以直接等价
    Canal_New_NUM=gg-1
    ALLOCATE(n_canal(Canal_NUM))
    j=0
    DO i=1,Element_NUM
        IF(Ele_Rel(i,2)==5)THEN
            j=j+1
            n_canal(j)=nnn(i)
        END IF
    END DO
    ALLOCATE(Ca_cl_New(Canal_New_NUM))
    k=0
    DO i=1,Canal_NUM
        DO j=1,n_canal(i)
            k=k+1
            Ca_cl_New(k)%Serial_NUM=k
            Ca_cl_New(k)%Name=Ca_cl(i)%Name
            Ca_cl_New(k)%BottomWidth=Ca_cl(i)%BottomWidth
            Ca_cl_New(k)%SideSlopes=Ca_cl(i)%SideSlopes
            Ca_cl_New(k)%Incoord=Ca_cl(i)%Incoord+(j-1)*(Ca_cl(i)%Outcoord-Ca_cl(i)%Incoord)/n_canal(i)
            Ca_cl_New(k)%Outcoord=Ca_cl(i)%Incoord+j*(Ca_cl(i)%Outcoord-Ca_cl(i)%Incoord)/n_canal(i)
            Ca_cl_New(k)%InInvertElev=Ca_cl(i)%InInvertElev+(j-1)*(Ca_cl(i)%OutInvertElev-Ca_cl(i)%InInvertElev)/n_canal(i)
            Ca_cl_New(k)%OutInvertElev=Ca_cl(i)%InInvertElev+j*(Ca_cl(i)%OutInvertElev-Ca_cl(i)%InInvertElev)/n_canal(i)
            Ca_cl_New(k)%InTopElev=Ca_cl(i)%InTopElev+(j-1)*(Ca_cl(i)%OutTopElev-Ca_cl(i)%InTopElev)/n_canal(i)
            Ca_cl_New(k)%OutTopElev=Ca_cl(i)%InTopElev+j*(Ca_cl(i)%OutTopElev-Ca_cl(i)%InTopElev)/n_canal(i)
            Ca_cl_New(k)%Roughness=Ca_cl(i)%Roughness
        END DO
    END DO

    DEALLOCATE(h)
    DEALLOCATE(Q)
    DEALLOCATE(Z)
    DEALLOCATE(vc)
    DEALLOCATE(L)
End 
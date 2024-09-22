!////渠道，元件代码 5
!////////分成两个部分，一个是恒定状态的，一个是非恒定状态的
!////初始系数
SUBROUTINE CanalInitialStageCalc(j_NUM,Waterdepth_flag,Serial,h2,Q2,Q1,h1,vvol)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_main
    !TYPE(Condition)::Condition_0
    !TYPE(Config)::Config_main
    REAL::Psi
    REAL::PsiR
    REAL::GRAV
    REAL::Waterdepth_flag
    
    !/////子程序与主调程序之间的数据交换
    INTEGER::j_NUM  !元件系统编号
    INTEGER::Serial !对应渠道的局部编号
    REAL::h2,Q2,Q1,h1   !分别表示循环过程中进口断面和出口断面的水深和流量
    INTEGER::nn  !元件分段个数
    INTEGER::j,i
    REAL::aa,ee,dhh_S
    REAL::vvol  !元件蓄量
    REAL::DeltaX    !空间步长

    !//////各断面水深、流量
    REAL,ALLOCATABLE::QQ_S(:)   !元件各断面流量
    REAL,ALLOCATABLE::hh_S(:)   !元件各断面水位

    !/////子程序内部计算涉及到的参数
    REAL::Len   !元件长度 !一般不用关键词做变量 wjb0415
    REAL::S0    !元件底坡

    !//////////存储各个小元件的参数，n个元件
    REAL,ALLOCATABLE::BottomWidth_s(:)    !各小元件底宽
    REAL,ALLOCATABLE::SideSlopes_s(:) !各小元件边坡系数
    REAL,ALLOCATABLE::Incoord_s(:) !各小元件入口里程
    REAL,ALLOCATABLE::Outcoord_s(:) !各小元件出口里程
    REAL,ALLOCATABLE::InInvertElev_s(:)   !各小元件入口底高程
    REAL,ALLOCATABLE::OutInvertElev_s(:)   !各小元件出口底高程
    REAL,ALLOCATABLE::InTopElev_s(:)   !各小元件入口顶高程
    REAL,ALLOCATABLE::OutTopElev_s(:)   !各小元件出口顶高程
    
    !/////////存储各断面的数据，n+1个断面
    REAL,ALLOCATABLE::BB(:) !各断面水面宽度
    REAL,ALLOCATABLE::Ar(:) !各断面过水面积
    REAL,ALLOCATABLE::Sf(:) !各断面摩擦坡度
    REAL,ALLOCATABLE::KK(:) !各断面流量模数
    REAL,ALLOCATABLE::L_wet(:)  !各断面湿周长度
    REAL,ALLOCATABLE::RR(:) !各断面水力半径
    !------------------
    INTEGER :: CroSec_TYPE
    !-----------------------------

    Psi=Config_main.Coeff_set.Psi
    PsiR=Config_main.Coeff_set.PsiR
    GRAV=Config_main.Coeff_set.GRAV
    
    DeltaX=Boundary_cal.Step*SQRT(GRAV*h2)
    !WRITE(*,*)DeltaX
    Len=(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord)*1000  

    Com_element_main.n(j_NUM)=int(Len/DeltaX)
    IF(Com_element_main.n(j_NUM)<1)THEN
        Com_element_main.n(j_NUM)=1
    END IF

    if (Com_element_main.Ca_cl(Serial).CroSecTYPE==2) Com_element_main.n(j_NUM)=1!#TU加，不规则断面的话不分段计算
    nn=Com_element_main.n(j_NUM)
    DeltaX=Len/nn    !空间步长
    ALLOCATE(BottomWidth_s(nn))  
    ALLOCATE(SideSlopes_s(nn))
    ALLOCATE(Incoord_s(nn))
    ALLOCATE(Outcoord_s(nn)) 
    ALLOCATE(InInvertElev_s(nn))   
    ALLOCATE(OutInvertElev_s(nn))  
    ALLOCATE(InTopElev_s(nn))   
    ALLOCATE(OutTopElev_s(nn)) 
    ALLOCATE(BB(nn+1))
    ALLOCATE(Ar(nn+1))
    ALLOCATE(Sf(nn+1))
    ALLOCATE(KK(nn+1))
    ALLOCATE(QQ_S(nn+1))
    ALLOCATE(hh_S(nn+1))
    ALLOCATE(L_wet(nn+1))
    ALLOCATE(RR(nn+1))

    DO j=1,nn
        BottomWidth_s(j)=Com_element_main.Ca_cl(Serial).BottomWidth
        SideSlopes_s(j)=Com_element_main.Ca_cl(Serial).SideSlopes
        Incoord_s(j)=Com_element_main.Ca_cl(Serial).Incoord+(j-1)*(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord)/nn
        Outcoord_s(j)=Com_element_main.Ca_cl(Serial).Incoord+j*(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord)/nn
        InInvertElev_s(j)=Com_element_main.Ca_cl(Serial).InInvertElev+(j-1)*(Com_element_main.Ca_cl(Serial).OutInvertElev-Com_element_main.Ca_cl(Serial).InInvertElev)/nn
        OutInvertElev_s(j)=Com_element_main.Ca_cl(Serial).InInvertElev+j*(Com_element_main.Ca_cl(Serial).OutInvertElev-Com_element_main.Ca_cl(Serial).InInvertElev)/nn
        InTopElev_s(j)=Com_element_main.Ca_cl(Serial).InTopElev+(j-1)*(Com_element_main.Ca_cl(Serial).OutTopElev-Com_element_main.Ca_cl(Serial).InTopElev)/nn
        OutTopElev_s(j)=Com_element_main.Ca_cl(Serial).InTopElev+j*(Com_element_main.Ca_cl(Serial).OutTopElev-Com_element_main.Ca_cl(Serial).InTopElev)/nn
    END DO
    
    S0=(Com_element_main.Ca_cl(Serial).InInvertElev-Com_element_main.Ca_cl(Serial).OutInvertElev)/((Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord)*1000)
    DO j=1,nn+1
        QQ_S(j)=Q1    
    END DO
    !--------------------wjb0414
    CroSec_TYPE=Com_element_main.Ca_cl(Serial).CroSecTYPE
    !-------------------------
    hh_S(nn+1)=h2
    hh_S(nn)=hh_S(nn+1)

   !断面计算，区分断面类型，wjb0414 增加一个断面类型参数和断面计算子函数
    DO j=nn,1,-1
        hh_S(j)=hh_S(j+1)
        DO i=1,200,1 !////200是最大循环数
            !j+1断面
            CALL Abr(Serial,CroSec_TYPE,hh_S(j+1),BottomWidth_s(j),SideSlopes_s(j),OutInvertElev_s(j),BB(j+1),Ar(j+1),L_wet(j+1),1,0)  !#TU
            !水力半径
            RR(j+1)=Ar(j+1)/L_wet(j+1)
            !流量模数K=CAsqrt(R)
            KK(j+1)=1/Com_element_main.Ca_cl(Serial).Roughness*Ar(j+1)**(5.0/3.0)*L_wet(j+1)**(-2.0/3.0) !同化时把糙率设置成相同的就可以了。
            Sf(j+1)=QQ_S(j+1)*ABS(QQ_S(j+1))/KK(j+1)**2.0
            !j断面
            CALL Abr(Serial,CroSec_TYPE,hh_S(j),BottomWidth_s(j),SideSlopes_s(j),inInvertElev_s(j),BB(j),Ar(j),L_wet(j),0,0)!#TU

            RR(j)=Ar(j)/L_wet(j)
            KK(j)=1/Com_element_main.Ca_cl(Serial).Roughness*Ar(j)**(5.0/3.0)*L_wet(j)**(-2.0/3.0)!同化时把糙率设置成相同的就可以了。
            Sf(j)=QQ_S(j)*ABS(QQ_S(j))/KK(j)**2.0

            !动量方程
            !Tu BottomWidth_s(j)+2*SideSlopes_s(j)*hh_S(j)这一整个改为BB(j)，避免不规则断面底宽为高程对的问题；边坡系数假设他影响不大，先不管
            aa=Psi/DeltaX*QQ_S(j)**2*BB(j)/(Ar(j)**3)-GRAV*Psi/DeltaX-2.0*Psi*(1.0-PsiR)*GRAV*Sf(j)/KK(j)*KK(j)/Ar(j)*(5.0/3.0*BB(j)-2.0/3.0*RR(j)*(2.0*SQRT(1.0+SideSlopes_s(j)**2)))!#TU改
            ee=-Psi/DeltaX/2.0*((QQ_S(j+1)/Ar(j+1))**2.0-(QQ_S(j)/Ar(j))**2.0)-(1.0-Psi)/DeltaX/2.0*((QQ_S(j+1)/Ar(j+1))**2.0-(QQ_S(j)/Ar(j))**2.0)-Psi*GRAV/DeltaX*(hh_S(j+1)-hh_S(j))-(1.0-Psi)*GRAV/DeltaX*(hh_S(j+1)-hh_S(j))-Psi*GRAV*(PsiR*Sf(j+1)+(1.0-PsiR)*Sf(j))-(1.0-Psi)*GRAV*(PsiR*Sf(j+1)+(1.0-PsiR)*Sf(j))+GRAV*S0

            !增量计算
            dhh_S=ee/aa
            hh_S(j)=hh_S(j)+dhh_S
            !TU问题出在这一步，一旦aa极接近0且小于0；
            !而ee大于0。那么hh_S是负的，后面一系列都错了

            !#TU加
            if (hh_S(j)<0) then
                write(*,*)"==============warning（恒定流迭代）======================"
                write(*,*) "在第",j_NUM,"个大元件（渠道）"
                write(*,*) "第",nn,"个小渠道上"
                write(*,*) "恒定流水深计算迭代过程中入口水位小于0，将其修改为3继续迭代"
                write(*,*) "如果此警告出现很多则代表迭代不过去，可以尝试修改3为其他值"
                write(*,*)"======================================================="
                hh_S(j)=3
            end if

            IF(ABS(dhh_S)<2.22*2.718281828459**(-16))THEN
                EXIT
            END IF
        END DO
        Waterdepth_flag=Waterdepth_flag+1
        Condition_0.h(Waterdepth_flag)=hh_S(j)
        Condition_0.Q(Waterdepth_flag)=QQ_S(j)
    END DO
    
    h1=hh_S(1)
    vvol=(h1+h2)/2.0*Com_element_main.Ca_cl(Serial).BottomWidth*1000.0*(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord) &
         &+Com_element_main.Ca_cl(Serial).SideSlopes*h1*h2*1000.0*(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord) &
         &+1.0/3.0*(Com_element_main.Ca_cl(Serial).SideSlopes*(h1-h2)**2)*1000.0*(Com_element_main.Ca_cl(Serial).Outcoord-Com_element_main.Ca_cl(Serial).Incoord)
    DEALLOCATE(BottomWidth_s)    
    DEALLOCATE(SideSlopes_s)
    DEALLOCATE(Incoord_s)
    DEALLOCATE(Outcoord_s) 
    DEALLOCATE(InInvertElev_s)   
    DEALLOCATE(OutInvertElev_s)  
    DEALLOCATE(InTopElev_s)   
    DEALLOCATE(OutTopElev_s)
    DEALLOCATE(BB)
    DEALLOCATE(Ar)
    DEALLOCATE(Sf)
    DEALLOCATE(KK)
    DEALLOCATE(QQ_S)
    DEALLOCATE(hh_S)
    DEALLOCATE(L_wet)
    DEALLOCATE(RR)
END

!////////////////非恒定流系数
Subroutine CanalDynamicCoeff(j_NUM,Serial,Qo2,Qo1,Zo2,Zo1,S2,T2,P2,V2,ele_vvol,Condition_cal)
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !TYPE(Boundary)::Boundary_cal
    !TYPE(Com_element)::Com_element_cal
    TYPE(Condition)::Condition_cal
    !TYPE(Config)::Config_main
    REAL::Theta
    REAL::GRAV
    REAL::RoughnessCoeff

    
    INTEGER::j_NUM  !元件系统编号    
    INTEGER::Serial !对应渠道的局部编号
    REAL::Qo2,Qo1,Zo2,Zo1 !分别表示n时刻j+1断面流量，j断面流量，j+1断面水位，j断面水位
    REAL::S2,T2,P2,V2   !追赶系数
    REAL::C,D,E,F,G,M   !线性方程组的系数
    REAL::Y1,Y2,Y3,Y4   !中间系数
    REAL::hh1,hh2   !水深
    REAL::BB1,BB2   !水面宽度
    REAL::Ar1,Ar2      !过水面积
    REAL::L_wet1,L_wet2     !湿周
    REAL::R1,R2 !水力半径
    REAL::C1,C2 !谢才系数
    REAL::u1,u2 !断面流速
    integer ::k             !#Tu计算时段，用于报错
    REAL::ele_vvol  !元件蓄量
    REAL::DeltaX    !空间步长
    
    Theta=Config_main.Coeff_set.Theta
    GRAV=Config_main.Coeff_set.GRAV
    RoughnessCoeff=Config_main.Coeff_set.RoughnessCoeff


    !/////需要用到的水力要素预处理
    hh1=Zo1-Com_element_cal.Ca_cl(Serial).InInvertElev
    hh2=Zo2-Com_element_cal.Ca_cl(Serial).OutInvertElev
    !#Tu-----------------------------------------以下均有修改
    if ((hh1<0).or.(hh2<0)) then
        write(*,*)
        write(*,*)
        write(*,*)"==============warning（非恒定流计算）==============="
        write(*,*) "在第",k,"时刻计算中"
        write(*,*) "第",j_NUM,"个计算元件（渠道）的"
        write(*,*) "非恒定流计算水深小于0"
        write(*,*) "说明上一时刻计算的水深不合理"
        write(*,*)"=================================================="
    end if

    CALL Abr(Serial,Com_element_cal.Ca_cl(Serial).CroSecTYPE,hh2,Com_element_cal.Ca_cl(Serial).BottomWidth,Com_element_cal.Ca_cl(Serial).SideSlopes,Com_element_cal.Ca_cl(Serial).OutInvertElev,BB2,Ar2,L_wet2,1,1)!TU#@底宽、边坡 如果是梯形断面，一整个渠道都是一样的，如果是不规则断面则用不到
    CALL Abr(Serial,Com_element_cal.Ca_cl(Serial).CroSecTYPE,hh1,Com_element_cal.Ca_cl(Serial).BottomWidth,Com_element_cal.Ca_cl(Serial).SideSlopes,Com_element_cal.Ca_cl(Serial).InInvertElev,BB1,Ar1,L_wet1,0,1)!TU#底高程则相反，只有不规则断面用的到。后三个是输出,最后两个是控制符

    if ((Ar1<0.001).or.(Ar2<0.001)) then
        write(*,*)
        write(*,*)
        write(*,*)"==============warning（非恒定流计算）==============="
        write(*,*) "在第",k,"时刻计算中"
        write(*,*) "第",j_NUM,"个计算元件（渠道）的"
        write(*,*) "面积为0"
        write(*,*) "说明有问题"
        write(*,*)"=================================================="
    end if

    !BB1=Com_element_cal.Ca_cl(Serial).BottomWidth+2.0*Com_element_cal.Ca_cl(Serial).SideSlopes*hh1
    !BB2=Com_element_cal.Ca_cl(Serial).BottomWidth+2.0*Com_element_cal.Ca_cl(Serial).SideSlopes*hh2
    DeltaX=(Com_element_cal.Ca_cl(Serial).Outcoord-Com_element_cal.Ca_cl(Serial).Incoord)*1000.0
    !Ar1=(BB1+Com_element_cal.Ca_cl(Serial).BottomWidth)*hh1/2.0
    !Ar2=(BB2+Com_element_cal.Ca_cl(Serial).BottomWidth)*hh2/2.0
    !L_wet1=Com_element_cal.Ca_cl(Serial).BottomWidth+2.0*hh1*SQRT(1+(Com_element_cal.Ca_cl(Serial).SideSlopes)**2.0)
    !L_wet2=Com_element_cal.Ca_cl(Serial).BottomWidth+2.0*hh2*SQRT(1+(Com_element_cal.Ca_cl(Serial).SideSlopes)**2.0)
    R1=Ar1/L_wet1
    R2=Ar2/L_wet2!#@TU
    u1=Qo1/Ar1
    u2=Qo2/Ar2
    !#Tu-----------------------------------------修改到此结束
    
    ele_vvol=(hh1+hh2)/2.0*Com_element_cal.Ca_cl(Serial).BottomWidth*1000.0*(Com_element_cal.Ca_cl(Serial).Outcoord-Com_element_cal.Ca_cl(Serial).Incoord) &
         +Com_element_cal.Ca_cl(Serial).SideSlopes*hh1*hh2*1000.0*(Com_element_cal.Ca_cl(Serial).Outcoord-Com_element_cal.Ca_cl(Serial).Incoord) &
         +1.0/3.0*(Com_element_cal.Ca_cl(Serial).SideSlopes*(hh1-hh2)**2)*1000.0*(Com_element_cal.Ca_cl(Serial).Outcoord-Com_element_cal.Ca_cl(Serial).Incoord)
    
    IF(Config_main.Assimilate_set.Control_sign==0)THEN     !模拟用
        C1=1.0/Com_element_cal.Ca_cl(Serial).Roughness*R1**(1.0/6.0)
        C2=1.0/Com_element_cal.Ca_cl(Serial).Roughness*R2**(1.0/6.0)    
    ELSE                        !同化用
        C1=1.0/RoughnessCoeff*R1**(1.0/6.0)
        C2=1.0/RoughnessCoeff*R2**(1.0/6.0)
    END IF

    !//////线性方程组系数赋值
    C=1.0/2.0*(BB1+BB2)*DeltaX/2.0/Boundary_cal.Step/Theta
    D=-(1.0-Theta)/Theta*(Qo2-Qo1)+C*(Zo2+Zo1)
    E=DeltaX/2.0/Theta/Boundary_cal.Step-u1+GRAV*ABS(u1)/2.0/Theta/C1**2/R1*DeltaX
    G=DeltaX/2.0/Theta/Boundary_cal.Step+u2+GRAV*ABS(u2)/2.0/Theta/C2**2/R2*DeltaX
    F=GRAV*1.0/2.0*(Ar1+Ar2)
    M=DeltaX/2.0/Theta/Boundary_cal.Step*(Qo2+Qo1)-(1-Theta)/Theta*(u2*Qo2-u1*Qo1)&
      -(1-Theta)/Theta*GRAV*1.0/2.0*(Ar1+Ar2)*(Zo2-Zo1)
    IF(Boundary_cal.Upsign==0)THEN
        !//////////中间系数
        Y1=Condition_cal.VV(j_NUM)+C
        Y2=F+E*Condition_cal.VV(j_NUM)
        Y3=D+Condition_cal.PP(j_NUM)
        Y4=M-E*Condition_cal.PP(j_NUM)

        !////////追赶系数
        S2=(G*Y3-Y4)/(Y1*G+Y2)
        T2=(C*G-F)/(Y1*G+Y2)
        P2=Y3-Y1*S2
        V2=C-Y1*T2
    ELSE
        !//////////中间系数
        Y1=D-C*Condition_cal.PP(j_NUM)
        Y2=M+F*Condition_cal.PP(j_NUM)
        Y3=1.0+C*Condition_cal.VV(j_NUM)
        Y4=E+F*Condition_cal.VV(j_NUM)

        !////////追赶系数
        S2=(C*Y2-F*Y1)/(F*Y3+C*Y4)
        T2=(C*G-F)/(F*Y3+C*Y4)
        P2=(Y1+Y3*S2)/C
        V2=(Y3*T2+1.0)/C
    END IF
END 
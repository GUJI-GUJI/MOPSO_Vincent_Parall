SUBROUTINE Steady
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    
    INTEGER::i,j,k
    !///////恒定流计算过程
    !//////给各断面水力要素赋初值
    k=0
    DO i=1,Element_NUM
        DO j=1,n(i)
            k=k+1
            QQ(k)=Q(i)
            hh(k)=h(i)+(j-1)*(h(i+1)-h(i))/n(i)
        END DO
    END DO
    QQ(c_element+1)=Q(Element_NUM+1)
    hh(c_element+1)=h(Element_NUM+1)
    !/////生成水动力模型恒定过程系数矩阵
    !///////上游边界处的系数
200 a(1)=0.0
    b(1)=0.0
    c(1)=0.0
    d(1)=1.0
    e(1)=Upstream(1)-QQ(1)
    !///////中间元件的系数
    DO j=1,c_element
       !/////先判断类型再调用对应的函数模块，生成相应的系数。
       IF(Ele_Rel_New(j,2)==5)THEN  !/////渠道元件
           Call CanalSteadyCoeff(Ele_Rel_New(j,3),hh(j),QQ(j),hh(j+1),QQ(j+1),a(2*j),b(2*j),c(2*j),d(2*j),e(2*j),a(2*j+1),b(2*j+1),c(2*j+1),d(2*j+1),e(2*j+1))
       ELSEIF(Ele_Rel_New(j,2)==7)THEN  !/////渐变段元件
           Call AreachangeSteadyCoeff(j,Ele_Rel_New(j,3),hh(j),QQ(j),hh(j+1),QQ(j+1),a(2*j),b(2*j),c(2*j),d(2*j),e(2*j),a(2*j+1),b(2*j+1),c(2*j+1),d(2*j+1),e(2*j+1))
       ELSEIF(Ele_Rel_New(j,2)==10)THEN !/////倒虹吸元件
           Call InvertedsiphonSteadyCoeff(j,Ele_Rel_New(j,3),hh(j),QQ(j),hh(j+1),QQ(j+1),a(2*j),b(2*j),c(2*j),d(2*j),e(2*j),a(2*j+1),b(2*j+1),c(2*j+1),d(2*j+1),e(2*j+1))
       ENDIF 
    END DO
    !////////下游边界处的系数
    a(Eq_NUM)=1.0
    b(Eq_NUM)=0.0
    c(Eq_NUM)=0.0
    d(Eq_NUM)=0.0
    e(Eq_NUM)=Downstream(1)-hh(c_element+1)
    
    !/////调用迭代算法求解线性方程组(追赶法)////////////
    !//////推求系数
    PP(1)=e(1)
    VV(1)=0.0
    DO j=1,c_element
        SS(j+1)=(e(2*j+1)*d(2*j+1)-e(2*j+1)*d(2*j)-(b(2*j)*d(2*j+1)-b(2*j+1)*d(2*j))*PP(j))/&
                (a(2*j)*d(2*j+1)-a(2*j+1)*d(2*j)-(b(2*j)*d(2*j+1)-b(2*j+1)*d(2*j))*VV(j))
        TT(j+1)=(c(2*j)*d(2*j+1)-c(2*j+1)*d(2*j))/&
                (a(2*j)*d(2*j+1)-a(2*j+1)*d(2*j)-(b(2*j)*d(2*j+1)-b(2*j+1)*d(2*j))*VV(j))
        PP(j+1)=((a(2*j+1)-b(2*j+1)*VV(j))*((e(2*j)-b(2*j)*PP(j))-(a(2*j)-b(2*j)*VV(j))*(e(2*j+1)-b(2*j+1)*PP(j)))/&
                ((a(2*j+1)-b(2*j+1)*VV(j))*d(2*j))-(a(2*j)-b(2*j)*VV(j))*d(2*j+1))
        VV(j+1)=((a(2*j+1)-b(2*j+1)*VV(j))*c(2*j)-(a(2*j)-b(2*j)*VV(j))*c(2*j+1))/&
                ((a(2*j+1)-b(2*j+1)*VV(j))*d(2*j)-(a(2*j)-b(2*j)*VV(j))*d(2*j+1))
    END DO
    !//////回代过程
    dh(c_element+1)=e(Eq_NUM)
    DO j=c_element,1,-1
        dh(j)=SS(j+1)-TT(j+1)*dh(j+1)
        dQ(j+1)=PP(j+1)-VV(j+1)*dh(j+1)
    END DO
    dQ(1)=e(1)
    !/////计算水深和流量
    DO j=1,c_element+1
        hh(j)=hh(j)+dh(j)
        QQ(j)=QQ(j)+dQ(j)
    END DO
    !////判断计算结果是否满足条件，如果不满足条件则重新进行计算
    DO j=1,c_element+1
        IF(ABS(dh(j))>0.0001.OR.ABS(dQ(j)>0.001))THEN
            GO TO 200
        END IF
    END DO
    
    !!////////存储本次计算结果，作为恒定流计算的初始条件
    ALLOCATE(hOld(c_element+1))
    ALLOCATE(QOld(c_element+1))
    ALLOCATE(ZOld(c_element+1))
    DO i=1,c_element+1
        hOld(i)=hh(i)
        QOld(i)=QQ(i)
    END DO
    
    !/////输出初始计算结果,按原有的计算元件输出/////////
    !//////计算各断面水位，即计算得到的水深加上各断面的底高程
    DO j=1,c_element,1
       IF(Ele_Rel_New(j,2)==5)THEN  !/////渠道元件
           ZOld(j)=Ca_cl_New(Ele_Rel_New(j,3))%InInvertElev+hOld(j)
           ZOld(j+1)=Ca_cl_New(Ele_Rel_New(j,3))%OutInvertElev+hOld(j+1)
       ELSEIF(Ele_Rel_New(j,2)==10)THEN !/////倒虹吸元件
           ZOld(j)=Insi_cl(Ele_Rel_New(j,3))%InInvertElev+hOld(j)
           ZOld(j+1)=Insi_cl(Ele_Rel_New(j,3))%OutInvertElev+hOld(j+1)
       ENDIF 
    END DO
    
    !/////输出恒定流计算水动力过程//////////
    OPEN(UNIT=21,FILE="SteadyRESULT.txt")
        WRITE(21,"(A7,3X,A7,3X,A7,3X,A7,3X)")"序号","水深","水位","流量"
        DO i=1,c_element+1
            WRITE(21,"(I7,3X,F7.3,3X,F7.3,3X,F7.3,3X)")i,hOld(i),ZOld(i),QOld(i)
        END DO
    CLOSE(21)
END
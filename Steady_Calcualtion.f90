SUBROUTINE Steady
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    
    INTEGER::i,j,k
    !///////�㶨���������
    !//////��������ˮ��Ҫ�ظ���ֵ
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
    !/////����ˮ����ģ�ͺ㶨����ϵ������
    !///////���α߽紦��ϵ��
200 a(1)=0.0
    b(1)=0.0
    c(1)=0.0
    d(1)=1.0
    e(1)=Upstream(1)-QQ(1)
    !///////�м�Ԫ����ϵ��
    DO j=1,c_element
       !/////���ж������ٵ��ö�Ӧ�ĺ���ģ�飬������Ӧ��ϵ����
       IF(Ele_Rel_New(j,2)==5)THEN  !/////����Ԫ��
           Call CanalSteadyCoeff(Ele_Rel_New(j,3),hh(j),QQ(j),hh(j+1),QQ(j+1),a(2*j),b(2*j),c(2*j),d(2*j),e(2*j),a(2*j+1),b(2*j+1),c(2*j+1),d(2*j+1),e(2*j+1))
       ELSEIF(Ele_Rel_New(j,2)==7)THEN  !/////�����Ԫ��
           Call AreachangeSteadyCoeff(j,Ele_Rel_New(j,3),hh(j),QQ(j),hh(j+1),QQ(j+1),a(2*j),b(2*j),c(2*j),d(2*j),e(2*j),a(2*j+1),b(2*j+1),c(2*j+1),d(2*j+1),e(2*j+1))
       ELSEIF(Ele_Rel_New(j,2)==10)THEN !/////������Ԫ��
           Call InvertedsiphonSteadyCoeff(j,Ele_Rel_New(j,3),hh(j),QQ(j),hh(j+1),QQ(j+1),a(2*j),b(2*j),c(2*j),d(2*j),e(2*j),a(2*j+1),b(2*j+1),c(2*j+1),d(2*j+1),e(2*j+1))
       ENDIF 
    END DO
    !////////���α߽紦��ϵ��
    a(Eq_NUM)=1.0
    b(Eq_NUM)=0.0
    c(Eq_NUM)=0.0
    d(Eq_NUM)=0.0
    e(Eq_NUM)=Downstream(1)-hh(c_element+1)
    
    !/////���õ����㷨������Է�����(׷�Ϸ�)////////////
    !//////����ϵ��
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
    !//////�ش�����
    dh(c_element+1)=e(Eq_NUM)
    DO j=c_element,1,-1
        dh(j)=SS(j+1)-TT(j+1)*dh(j+1)
        dQ(j+1)=PP(j+1)-VV(j+1)*dh(j+1)
    END DO
    dQ(1)=e(1)
    !/////����ˮ�������
    DO j=1,c_element+1
        hh(j)=hh(j)+dh(j)
        QQ(j)=QQ(j)+dQ(j)
    END DO
    !////�жϼ������Ƿ�����������������������������½��м���
    DO j=1,c_element+1
        IF(ABS(dh(j))>0.0001.OR.ABS(dQ(j)>0.001))THEN
            GO TO 200
        END IF
    END DO
    
    !!////////�洢���μ���������Ϊ�㶨������ĳ�ʼ����
    ALLOCATE(hOld(c_element+1))
    ALLOCATE(QOld(c_element+1))
    ALLOCATE(ZOld(c_element+1))
    DO i=1,c_element+1
        hOld(i)=hh(i)
        QOld(i)=QQ(i)
    END DO
    
    !/////�����ʼ������,��ԭ�еļ���Ԫ�����/////////
    !//////���������ˮλ��������õ���ˮ����ϸ�����ĵ׸߳�
    DO j=1,c_element,1
       IF(Ele_Rel_New(j,2)==5)THEN  !/////����Ԫ��
           ZOld(j)=Ca_cl_New(Ele_Rel_New(j,3))%InInvertElev+hOld(j)
           ZOld(j+1)=Ca_cl_New(Ele_Rel_New(j,3))%OutInvertElev+hOld(j+1)
       ELSEIF(Ele_Rel_New(j,2)==10)THEN !/////������Ԫ��
           ZOld(j)=Insi_cl(Ele_Rel_New(j,3))%InInvertElev+hOld(j)
           ZOld(j+1)=Insi_cl(Ele_Rel_New(j,3))%OutInvertElev+hOld(j+1)
       ENDIF 
    END DO
    
    !/////����㶨������ˮ��������//////////
    OPEN(UNIT=21,FILE="SteadyRESULT.txt")
        WRITE(21,"(A7,3X,A7,3X,A7,3X,A7,3X)")"���","ˮ��","ˮλ","����"
        DO i=1,c_element+1
            WRITE(21,"(I7,3X,F7.3,3X,F7.3,3X,F7.3,3X)")i,hOld(i),ZOld(i),QOld(i)
        END DO
    CLOSE(21)
END
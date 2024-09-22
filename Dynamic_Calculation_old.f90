    SUBROUTINE Dynamic_Cal(C_Opendegree_FF,ST_STATE_FF,DS_STATE_FF,Pool_Volume,QControlGate_FF)!
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    INTEGER::i,j,k,mm,kk,LLL
    real::c11,c22,c21,hh
    REAL::zzaverage,QQaverage
    INTEGER::si_n   !�Ż��洢ʱ��ѭ������
    REAL::random    !�����
    DOUBLE PRECISION r0
    Integer::pos_sg,DM_sg !��¼����բλ��
    INTEGER::ele_sg    !��¼��Ԫ����λ�ã������Ŷգ�
    Integer::st,sl  !�洢����բ��Ϣ���
    !/////////////�Ǻ㶨������
    !δ֪���洢
    REAL,ALLOCATABLE::ZZ(:,:)
    REAL,ALLOCATABLE::QQ(:,:)
    !��һ����һֱ���洢
    REAL,ALLOCATABLE::ZZOLD(:,:)
    REAL,ALLOCATABLE::QQOLD(:,:)
    !׷��ϵ��
    REAL,ALLOCATABLE::SS(:)
    REAL,ALLOCATABLE::TT(:)
    REAL,ALLOCATABLE::PP(:)
    REAL,ALLOCATABLE::VV(:)
    !��̬������ˮ�����飨�£�
    REAL,ALLOCATABLE::h_new(:)
    INTEGER,ALLOCATABLE::SG_POSTION(:),DM_POSTION(:),ELE_POSTION(:)
    REAL::C_Opendegree_FF(Controlgate_NUM,Calculate_NUM)
    REAL::QControlGate_FF((End_Time-Start_Time)/3600+1,2*Controlgate_NUM)
    REAL::ST_STATE_FF((End_Time-Start_Time)/3600+1,2*Controlgate_NUM)
    REAL::DS_STATE_FF((End_Time-Start_Time)/3600+1)
    REAL::Pool_Volume(Controlgate_NUM+1,Calculate_NUM)
    REAL,ALLOCATABLE::ElE_Volume(:,:)
    REAL,ALLOCATABLE::Total_Volume(:)
    !    ALLOCATE(C_Opendegree_FF(Controlgate_NUM,Calculate_NUM))
    !    ALLOCATE(ST_STATE_FF((End_Time-Start_Time)/3600+1,2*Controlgate_NUM))
    !    ALLOCATE(DS_STATE_FF((End_Time-Start_Time)/3600+1))
    !//////////////////////////////////
    ALLOCATE(h_new(c_element+1))
    ALLOCATE(ZZ(c_element+1,NUMMember))
    ALLOCATE(QQ(c_element+1,NUMMember))
    ALLOCATE(ZZOLD(c_element+1,NUMMember))
    ALLOCATE(QQOLD(c_element+1,NUMMember))
    ALLOCATE(SG_POSTION(Controlgate_NUM))
    ALLOCATE(DM_POSTION(Controlgate_NUM+2))
    ALLOCATE(ELE_POSTION(Element_NUM-Bridge_NUM))
    ALLOCATE(Dynamic_Result_Zup(Controlgate_NUM,Calculate_NUM))
    ALLOCATE(Dynamic_Result_Zdown(Controlgate_NUM,Calculate_NUM))
    ALLOCATE(Dynamic_Result_Q(Controlgate_NUM,Calculate_NUM))
    ALLOCATE(Dynamic_Result_Zbd(1,Calculate_NUM))
    ALLOCATE(Dynamic_Result_Qbd(1,Calculate_NUM))
    ALLOCATE(Dynamic_Result_ZZ(Element_NUM-Bridge_NUM,Calculate_NUM))
    ALLOCATE(Dynamic_Result_QQ(Element_NUM-Bridge_NUM,Calculate_NUM))
    ALLOCATE(ElE_Volume(c_element,NUMMember))
    ALLOCATE(Total_Volume(Calculate_NUM))
    !///////ȷ��բ�����ڵ�λ��//////
    pos_sg=1
    DO j=1,c_element
        IF(Ele_Rel_New(j,2)==6)THEN
            SG_POSTION(pos_sg)=j
            pos_sg=pos_sg+1
        ENDIF
    END DO

    !///////ȷ��բ��+��λ�������ڵ�λ��//////
    DM_sg=2
    DM_POSTION(1)=1
    DO j=2,c_element
        IF(Ele_Rel_New(j,2)==6)THEN
            DM_POSTION(DM_sg)=j
            DM_sg=DM_sg+1
        END IF
    END DO
    DM_POSTION(Controlgate_NUM+2)=c_element

    !///////ȷ�����������ڵ�λ��//////
    ele_sg=2
    ELE_POSTION(1)=1
    DO j=2,c_element
        IF(Ele_Rel_New(j,2)==4.OR.Ele_Rel_New(j,2)==6.OR.Ele_Rel_New(j,2)==7.OR.Ele_Rel_New(j,2)==10.OR.Ele_Rel_New(j,2)==12)THEN
            ELE_POSTION(ele_sg)=j
            ele_sg=ele_sg+1
        ELSEIF(Ele_Rel_New(j-1,2)/=5.AND.Ele_Rel_New(j,2)==5)THEN
            ELE_POSTION(ele_sg)=j
            ele_sg=ele_sg+1
        ENDIF
    END DO


    DO j=c_element+1,1,-1
        h_new(j)=Waterdepth(c_element+2-j,1)
        QQ(j,1)=Waterdepth(c_element+2-j,2)
    END DO

    DO mm=1,NUMMember
        DO j=1,c_element,1
            IF(Ele_Rel_New(j,2)==5)THEN  !/////����Ԫ��
                ZZ(j,mm)=Ca_cl_New(Ele_Rel_New(j,3))%InInvertElev+h_new(j)
                ZZ(j+1,mm)=Ca_cl_New(Ele_Rel_New(j,3))%OutInvertElev+h_new(j+1)
            ELSEIF(Ele_Rel_New(j,2)==10)THEN !/////������Ԫ��
                ZZ(j,mm)=Insi_cl(Ele_Rel_New(j,3))%InInvertElev+h_new(j)
                ZZ(j+1,mm)=Insi_cl(Ele_Rel_New(j,3))%OutInvertElev+h_new(j+1)
            ENDIF
        END DO
    END DO
    IF(Boundary_sign_UP==1)THEN
        ZZ(1,1)=Upstream(1)
    ENDIF

    !////////�ܵļ������Ϊ��Calculate_NUM,�������ļ��м������
    si_n=1
    st=0
    DO k=1,Calculate_NUM
        DO mm=1,NUMMember
            DO i=1,c_element+1,1
                QQOLD(i,mm)=QQ(i,mm)
                ZZOLD(i, mm)=ZZ(i,mm)
            END DO
            !׷��ϵ��
            ALLOCATE(SS(c_element+1))
            ALLOCATE(TT(c_element+1))
            ALLOCATE(PP(c_element+1))
            ALLOCATE(VV(c_element+1))


            !***************����ϵ��************
            IF(Boundary_sign_UP==0)THEN        !���α߽������ֱ�ѡ��ˮλ������.0��ʾѡ�������߽磬1��ʾ����Ϊˮλ�߽�
                !/////////���α߽�ϵ��
                PP(1)=Upstream(k)
                QQ(1,mm)=Upstream(k)
                VV(1)=0.0
                !/////////�м�Ԫ��ϵ��
                DO j=1,c_element
                    !/////���ж������ٵ��ö�Ӧ�ĺ���ģ�飬������Ӧ��ϵ����
                    IF(Ele_Rel_New(j,2)==5)THEN  !/////����Ԫ��
                        Call CanalDynamicCoeff(j,Ele_Rel_New(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==4)THEN     !////��ˮ��Ԫ��
                        Call SideOutFlowDynamicCoeff(j,Ele_Rel_New(j,3),C_OutDischarge(Ele_Rel_New(j,3),k),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==6)THEN     !////����բԪ��
                        Call ControlGateDynamicCoeff(j,Ele_Rel_New(j,3),C_Opendegree_FF(Ele_Rel_New(j,3),k),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==7)THEN  !/////�����Ԫ��
                        Call AreachangeDynamicCoeff(j,Ele_Rel_New(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==10)THEN !/////������Ԫ��
                        Call InvertedsiphonDynamicCoeff(j,Ele_Rel_New(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==11)THEN !/////�Ŷ�Ԫ��
                        Call BridgeDynamicCoeff(j,Ele_Rel_New(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==12)THEN    !/////��վԪ��
                        Call PumpingstationDynamicCoeff(j,Ele_Rel_New(j,3),C_Bladeangle_cl(Ele_Rel_New(j,3))%Each_PS(:,k),PumpST_cl(Ele_Rel_New(j,3))%UNIT_NUM,QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ENDIF
                END DO
                !///////////���α߽�����
                IF(Boundary_sign_DOWN==1)THEN   !///////////ˮλ�߽�1
                    !ZZ(c_element+1,mm)=Downstream(k)+Ca_cl_New(Ele_Rel_New(c_element,3))%OutInvertElev
                    ZZ(c_element+1,mm)=Downstream(k)
                ELSEIF(Boundary_sign_DOWN==2)THEN   !///////////ˮλ~������ϵ2
                    !m_coeff=-0.001*(Downstream(k)/1000)**2-0.0673*(Downstream(k)/1000)+0.4323
                    m_coeff=0.35     !ȡ�����ɳ�������ֵ
                    PPQ=m_coeff*21.0*Downstream(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Ca_cl(Canal_NUM)%OutInvertElev))-1.0/(2.0*m_coeff*21.0*Downstream(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Ca_cl(Canal_NUM)%OutInvertElev))*(ZZOLD(c_element+1,mm)-Ca_cl(Canal_NUM)%OutInvertElev))
                    VVQ=-1.0/(2.0*m_coeff*21.0*Downstream(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Ca_cl(Canal_NUM)%OutInvertElev)))
                    ZZ(c_element+1,mm)=(PP(c_element+1)-PPQ)/(VV(c_element+1)-VVQ)
                ENDIF

                !*******�ش����ˮ��������********��
                DO j=c_element,1,-1
                    ZZ(j,mm)=SS(j+1)-TT(j+1)*ZZ(j+1,mm)
                    QQ(j+1,mm)=PP(j+1)-VV(j+1)*ZZ(j+1,mm)
                END DO
            ELSE    !����Ϊˮλ�߽�Boundary_sign_UP==1
                !/////////���α߽�ϵ��
                PP(1)=Upstream(k)
                ZZ(1,mm)=Upstream(k)
                VV(1)=0.0
                !/////////�м�Ԫ��ϵ��
                DO j=1,c_element
                    !/////���ж������ٵ��ö�Ӧ�ĺ���ģ�飬������Ӧ��ϵ����
                    IF(Ele_Rel_New(j,2)==5)THEN  !/////����Ԫ��
                        Call CanalDynamicCoeff(j,Ele_Rel_New(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==4)THEN     !////��ˮ��Ԫ��
                        Call SideOutFlowDynamicCoeff(j,Ele_Rel_New(j,3),C_OutDischarge(Ele_Rel_New(j,3),k),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==6)THEN     !////����բԪ��
                        Call ControlGateDynamicCoeff(j,Ele_Rel_New(j,3),C_Opendegree_FF(Ele_Rel_New(j,3),k),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==7)THEN  !/////�����Ԫ��
                        Call AreachangeDynamicCoeff(j,Ele_Rel_New(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==10)THEN !/////������Ԫ��
                        Call InvertedsiphonDynamicCoeff(j,Ele_Rel_New(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==11)THEN !/////�Ŷ�Ԫ��
                        Call BridgeDynamicCoeff(j,Ele_Rel_New(j,3),QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ELSEIF(Ele_Rel_New(j,2)==12)THEN    !/////��վԪ��
                        Call PumpingstationDynamicCoeff(j,Ele_Rel_New(j,3),C_Bladeangle_cl(Ele_Rel_New(j,3))%Each_PS(:,k),PumpST_cl(Ele_Rel_New(j,3))%UNIT_NUM,QQOLD(j+1,mm),QQOLD(j,mm),ZZOLD(j+1,mm),ZZOLD(j,mm),SS(j),TT(j),PP(j),VV(j),SS(j+1),TT(j+1),PP(j+1),VV(j+1),ElE_Volume(j,mm))
                    ENDIF
                END DO
                !///////////���α߽�����
                IF(Boundary_sign_DOWN==0)THEN   !///////////�����߽�0
                    QQ(c_element+1,mm)=Downstream(k)
                ELSEIF(Boundary_sign_DOWN==1)THEN   !///////////ˮλ�߽�1
                    QQ(c_element+1,mm)=(PP(c_element+1)-Downstream(k))/VV(c_element+1)
                ELSEIF(Boundary_sign_DOWN==2)THEN   !///////////ˮλ~������ϵ�߽�2
                    !m_coeff=-0.001*(Downstream(k)/1000)**2-0.0673*(Downstream(k)/1000)+0.4323
                    m_coeff=0.35      !ȡ�����ɳ�������ֵ
                    PPQ=m_coeff*21.0*Downstream(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Ca_cl(Canal_NUM)%OutInvertElev))-1.0/(2.0*m_coeff*21.0*Downstream(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Ca_cl(Canal_NUM)%OutInvertElev))*(ZZOLD(c_element+1,mm)-Ca_cl(Canal_NUM)%OutInvertElev))
                    VVQ=-1.0/(2.0*m_coeff*21.0*Downstream(k)/1000*SQRT(2*9.81*(ZZOLD(c_element+1,mm)-Ca_cl(Canal_NUM)%OutInvertElev)))
                    QQ(c_element+1,mm)=(PPQ-VVQ*PP(c_element+1))/(1-VVQ*VV(c_element+1))
                ENDIF
                !*******�ش����ˮ��������********��
                DO j=c_element,1,-1
                    QQ(j,mm)=SS(j+1)-TT(j+1)*QQ(j+1,mm)
                    ZZ(j+1,mm)=PP(j+1)-VV(j+1)*QQ(j+1,mm)
                END DO

            END IF



            IF(mod((k-1)*DeltaT,3600)==0)THEN
                st=st+1
                !���α߽紦��ˮλ����
                DS_STATE_FF(st)=ZZ(c_element+1,mm)
                DO DM_sg=1,Controlgate_NUM
                    ST_STATE_FF(st,2*DM_sg-1)=ZZ(SG_POSTION(DM_sg),mm)
                    ST_STATE_FF(st,2*DM_sg)=ZZ(SG_POSTION(DM_sg)+1,mm)
                    QControlGate_FF(st,2*DM_sg)=QQ(SG_POSTION(DM_sg)+1,mm)
                ENDDO
            END IF
            !***���������ܺ�*******
            DO i=1,Controlgate_NUM+1
                Pool_Volume(i,k)=0.0
                DO j=DM_POSTION(i),DM_POSTION(i+1)-1
                    Pool_Volume(i,k)=Pool_Volume(i,k)+ElE_Volume(j,mm)
                END DO
            END DO
            !Pool_Volume(Controlgate_NUM+1,k)=Pool_Volume(Controlgate_NUM+1,k)+ElE_Volume(c_element,mm)

            !***�����ܺ�*******
            Total_Volume(k)=0.0
            DO i=1,Controlgate_NUM+1
                Total_Volume(k)=Total_Volume(k)+Pool_volume(i,k)
            END DO


            !*******��մ洢�ռ�WRITE(24,*)
            DEALLOCATE(SS)
            DEALLOCATE(TT)
            DEALLOCATE(PP)
            DEALLOCATE(VV)
        END DO
    END DO

    !////////////////////////////////////////////////////////////////////////////////////////////////

    !///////�ͷſռ�
    DEALLOCATE(h_new)
    !DEALLOCATE(ZZ)
    !DEALLOCATE(QQ)
    DEALLOCATE(ZZOLD)
    DEALLOCATE(QQOLD)
    DEALLOCATE(SG_POSTION)
    END
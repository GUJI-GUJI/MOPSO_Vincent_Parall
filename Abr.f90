SUBROUTINE ABR(SEC_ID,SEC_TYPE,ZH,b0,Side_Slope,ZDMIN,BW,A,X,inout_flag,dynamic_flag)   !#tu�����������Ϊ �����š����ͣ�ˮ��׿����£�������������Ϊˮ��������ʪ��,�����ڶ����ʶ��tu�ӣ�0Ϊ������ڣ�1Ϊ�������ڣ����Ǻ㶨�������ʶ
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !--------------------------------
    !TYPE(Com_element)::Com_element_main
    INTEGER,INTENT(IN) :: SEC_ID,SEC_TYPE
    REAL,INTENT(IN) :: ZH,b0,Side_Slope,ZDMIN  !----- ZDMIN������͸߳�
    REAL,INTENT(OUT) :: BW, A,X
    INTEGER :: ND
    REAL :: ZD1,ZD2,BD1,BD2,BD,H1,H2,ZDI
    INTEGER :: JJJ
    REAL :: ZMAX !,ZMIN   !�����������ˮλ�����ˮλ
    REAL :: ZW    !ˮλ����ˮ��߳�
    integer :: new_SEC_ID!#TU��new_SEC_ID��SEC_ID+1,ʹ����ȡ��GCDʱ��ZDMIN��Ӧ����
    Integer :: inout_flag!#TU,0Ϊ������ڣ�1Ϊ��������
    Integer :: dynamic_flag!#TU��Ϊ0�Ǻ㶨�����㣬��Ca_cl��������Ca_cl_new

    new_SEC_ID=SEC_ID+1;!#tu
    SELECTCASE (SEC_TYPE)
    CASE(0)   !����
        BW=b0+2*Side_Slope*ZH
        A=(BW+b0)*ZH/2.0
        X=b0+2*ZH*SQRT(1+Side_Slope**2.0)
    CASE(1)   !Բ�ζ���
        IF(ZH>=2*b0) THEN    !Բ�ζ���b0Ϊ�뾶
            WRITE(*,"('Բ������')")
            BW=0.0
            A=PI*b0**2.0
            X=PI*b0
        ELSEIF(ZH>=b0) THEN
            BW=2*SQRT(b0**2.0-(ZH-b0)**2.0)
            A=PI*b0**2.0*(1-ACOS((ZH-b0)/b0)/PI)+0.5*BW*(ZH-b0)
            X=2*PI*b0*(1-ACOS((ZH-b0)/b0)/PI)
        ELSE
            BW=2*SQRT(b0**2.0-(b0-ZH)**2.0)
            A=PI*b0**2.0*ACOS((b0-ZH)/b0)/PI-0.5*BW*(b0-ZH)
            X=2*PI*b0*ACOS((b0-ZH)/b0)
        END IF
    CASE(2) !���������, ���������ʱ����ˮλ
        !------------------------------------------#tu����һ���ζ����޸ģ�����ʵ�ֱַ��ȡ�����ڵĶ��������Ʋ���������
        A=0.0
        BW=0.0
        X=0.0
        ZMAX=-1000.0
        ZW=ZH+ZDMIN
        if (dynamic_flag==0) then!�㶨������
            if (inout_flag==0) then!��ȡ��ڶ���
                ND=Com_element_main.Ca_cl(SEC_ID).CDS_in
                !-----------------------ȷ��������ͺ���ߵ�
                DO JJJ=1,ND!Ѱ�Ҹ̵߳��е����ֵ
                    IF(Com_element_main.Ca_cl(SEC_ID).GCD_in(2,JJJ)>ZMAX) ZMAX=Com_element_main.Ca_cl(SEC_ID).GCD_in(2,JJJ)
                END DO
                !---------------
                DO JJJ=1,ND-1
                    ZD1=Com_element_main.Ca_cl(SEC_ID).GCD_in(2,JJJ)
                    BD1=Com_element_main.Ca_cl(SEC_ID).GCD_in(1,JJJ)
                    ZD2=Com_element_main.Ca_cl(SEC_ID).GCD_in(2,JJJ+1)
                    BD2=Com_element_main.Ca_cl(SEC_ID).GCD_in(1,JJJ+1)
                    ZDI=MIN(ZD1,ZD2)
                    !--------------
                    !��������ˮλ���������ĸ߳̾��ͣ���ѭ������J��1�����ִ��ѭ��
                    IF(ZW<=ZDI)THEN
                        CYCLE	 !CYCLE�������ֹ����ѭ����ʹѭ��������1�ټ���ִ��ѭ��
                    ELSE
                        H1=ZW-ZD1
                        H2=ZW-ZD2
                        BD=BD2-BD1
                        IF(H1<0)THEN
                            H1=0
                            BD=H2/(ZD1-ZD2)*BD
                        ELSEIF(H2<0)THEN
                            H2=0
                            BD=H1/(ZD2-ZD1)*BD
                        ENDIF
                        A=A+0.5*(H1+H2)*BD
                        X=X+(ABS(H1-H2)**2+BD**2)**0.5
                        BW=BW+BD
                    ENDIF
                END DO !JJJ

            else !TU��ȡ���ڶ���
                ND=Com_element_main.Ca_cl(SEC_ID).CDS_out
                !-----------------------ȷ��������ͺ���ߵ�
                DO JJJ=1,ND!Ѱ�Ҹ̵߳��е����ֵ
                    IF(Com_element_main.Ca_cl(SEC_ID).GCD_out(2,JJJ)>ZMAX) ZMAX=Com_element_main.Ca_cl(SEC_ID).GCD_out(2,JJJ)
                END DO
                !---------------
                DO JJJ=1,ND-1
                    ZD1=Com_element_main.Ca_cl(SEC_ID).GCD_out(2,JJJ)
                    BD1=Com_element_main.Ca_cl(SEC_ID).GCD_out(1,JJJ)
                    ZD2=Com_element_main.Ca_cl(SEC_ID).GCD_out(2,JJJ+1)
                    BD2=Com_element_main.Ca_cl(SEC_ID).GCD_out(1,JJJ+1)
                    ZDI=MIN(ZD1,ZD2)
                    !--------------
                    !��������ˮλ���������ĸ߳̾��ͣ���ѭ������J��1�����ִ��ѭ��
                    IF(ZW<=ZDI)THEN
                        CYCLE	 !CYCLE�������ֹ����ѭ����ʹѭ��������1�ټ���ִ��ѭ��
                    ELSE
                        H1=ZW-ZD1
                        H2=ZW-ZD2
                        BD=BD2-BD1
                        IF(H1<0)THEN
                            H1=0
                            BD=H2/(ZD1-ZD2)*BD
                        ELSE IF(H2<0)THEN
                            H2=0
                            BD=H1/(ZD2-ZD1)*BD
                        END IF
                        A=A+0.5*(H1+H2)*BD
                        X=X+(ABS(H1-H2)**2+BD**2)**0.5
                        BW=BW+BD
                    END IF
                END DO !JJJ
            END IF
        else!�Ǻ㶨������
            if (inout_flag==0) then!��ȡ��ڶ���
                ND=Com_element_main.Ca_cl(SEC_ID).CDS_in
                !-----------------------ȷ��������ͺ���ߵ�
                DO JJJ=1,ND!Ѱ�Ҹ̵߳��е����ֵ
                    IF(Com_element_main.Ca_cl(SEC_ID).GCD_in(2,JJJ)>ZMAX) ZMAX=Com_element_main.Ca_cl(SEC_ID).GCD_in(2,JJJ)
                END DO
                !---------------
                DO JJJ=1,ND-1
                    ZD1=Com_element_main.Ca_cl(SEC_ID).GCD_in(2,JJJ)
                    BD1=Com_element_main.Ca_cl(SEC_ID).GCD_in(1,JJJ)
                    ZD2=Com_element_main.Ca_cl(SEC_ID).GCD_in(2,JJJ+1)
                    BD2=Com_element_main.Ca_cl(SEC_ID).GCD_in(1,JJJ+1)
                    ZDI=MIN(ZD1,ZD2)
                    !--------------
                    !��������ˮλ���������ĸ߳̾��ͣ���ѭ������J��1�����ִ��ѭ��
                    IF(ZW<=ZDI)THEN
                        CYCLE	 !CYCLE�������ֹ����ѭ����ʹѭ��������1�ټ���ִ��ѭ��
                    ELSE
                        H1=ZW-ZD1
                        H2=ZW-ZD2
                        BD=BD2-BD1
                        IF(H1<0)THEN
                            H1=0
                            BD=H2/(ZD1-ZD2)*BD
                        ELSEIF(H2<0)THEN
                            H2=0
                            BD=H1/(ZD2-ZD1)*BD
                        ENDIF
                        A=A+0.5*(H1+H2)*BD
                        X=X+(ABS(H1-H2)**2+BD**2)**0.5
                        BW=BW+BD
                    ENDIF
                END DO !JJJ
            else !TU��ȡ���ڶ���
                ND=Com_element_main.Ca_cl(SEC_ID).CDS_out
                !-----------------------ȷ��������ͺ���ߵ�
                DO JJJ=1,ND!Ѱ�Ҹ̵߳��е����ֵ
                    IF(Com_element_main.Ca_cl(SEC_ID).GCD_out(2,JJJ)>ZMAX) ZMAX=Com_element_main.Ca_cl(SEC_ID).GCD_out(2,JJJ)
                END DO
                !---------------
                DO JJJ=1,ND-1
                    ZD1=Com_element_main.Ca_cl(SEC_ID).GCD_out(2,JJJ)
                    BD1=Com_element_main.Ca_cl(SEC_ID).GCD_out(1,JJJ)
                    ZD2=Com_element_main.Ca_cl(SEC_ID).GCD_out(2,JJJ+1)
                    BD2=Com_element_main.Ca_cl(SEC_ID).GCD_out(1,JJJ+1)
                    ZDI=MIN(ZD1,ZD2)
                    !--------------
                    !��������ˮλ���������ĸ߳̾��ͣ���ѭ������J��1�����ִ��ѭ��
                    IF(ZW<=ZDI)THEN
                        CYCLE	 !CYCLE�������ֹ����ѭ����ʹѭ��������1�ټ���ִ��ѭ��
                    ELSE
                        H1=ZW-ZD1
                        H2=ZW-ZD2
                        BD=BD2-BD1
                        IF(H1<0)THEN
                            H1=0
                            BD=H2/(ZD1-ZD2)*BD
                        ELSE IF(H2<0)THEN
                            H2=0
                            BD=H1/(ZD2-ZD1)*BD
                        END IF
                        A=A+0.5*(H1+H2)*BD
                        X=X+(ABS(H1-H2)**2+BD**2)**0.5
                        BW=BW+BD
                    END IF
                END DO !JJJ
            END IF
        end if
        !------------------------------------------#Tu�޸ĵ���Ϊֹ
        !------------------------------------------
        CASE DEFAULT
    END SELECT
    !--------------------------------------------------
    RETURN
    END SUBROUTINE
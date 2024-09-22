SUBROUTINE ABR(SEC_ID,SEC_TYPE,ZH,b0,Side_Slope,ZDMIN,BW,A,X,inout_flag,dynamic_flag)   !#tu传入参数依次为 断面编号、类型，水深，底宽，边坡，传出参数依次为水面宽，面积，湿周,进出口断面标识（tu加，0为渠道入口，1为渠道出口），非恒定流计算标识
    USE GLOBAL_VARIABLE
    
    IMPLICIT NONE
    !--------------------------------
    !TYPE(Com_element)::Com_element_main
    INTEGER,INTENT(IN) :: SEC_ID,SEC_TYPE
    REAL,INTENT(IN) :: ZH,b0,Side_Slope,ZDMIN  !----- ZDMIN断面最低高程
    REAL,INTENT(OUT) :: BW, A,X
    INTEGER :: ND
    REAL :: ZD1,ZD2,BD1,BD2,BD,H1,H2,ZDI
    INTEGER :: JJJ
    REAL :: ZMAX !,ZMIN   !断面允许最高水位、最低水位
    REAL :: ZW    !水位，即水面高程
    integer :: new_SEC_ID!#TU，new_SEC_ID是SEC_ID+1,使得在取用GCD时与ZDMIN对应断面
    Integer :: inout_flag!#TU,0为渠道入口，1为渠道出口
    Integer :: dynamic_flag!#TU，为0是恒定流计算，用Ca_cl；否则用Ca_cl_new

    new_SEC_ID=SEC_ID+1;!#tu
    SELECTCASE (SEC_TYPE)
    CASE(0)   !梯形
        BW=b0+2*Side_Slope*ZH
        A=(BW+b0)*ZH/2.0
        X=b0+2*ZH*SQRT(1+Side_Slope**2.0)
    CASE(1)   !圆形断面
        IF(ZH>=2*b0) THEN    !圆形断面b0为半径
            WRITE(*,"('圆管满流')")
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
    CASE(2) !不规则断面, 不规则断面时输入水位
        !------------------------------------------#tu以下一整段都有修改，尝试实现分别读取进出口的断面来控制不规则渠道
        A=0.0
        BW=0.0
        X=0.0
        ZMAX=-1000.0
        ZW=ZH+ZDMIN
        if (dynamic_flag==0) then!恒定流计算
            if (inout_flag==0) then!读取入口断面
                ND=Com_element_main.Ca_cl(SEC_ID).CDS_in
                !-----------------------确定断面最低和最高点
                DO JJJ=1,ND!寻找高程点中的最大值
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
                    !若给定的水位比两个测点的高程均低，则循环变量J加1后继续执行循环
                    IF(ZW<=ZDI)THEN
                        CYCLE	 !CYCLE语句是终止本次循环，使循环变量加1再继续执行循环
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

            else !TU读取出口断面
                ND=Com_element_main.Ca_cl(SEC_ID).CDS_out
                !-----------------------确定断面最低和最高点
                DO JJJ=1,ND!寻找高程点中的最大值
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
                    !若给定的水位比两个测点的高程均低，则循环变量J加1后继续执行循环
                    IF(ZW<=ZDI)THEN
                        CYCLE	 !CYCLE语句是终止本次循环，使循环变量加1再继续执行循环
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
        else!非恒定流计算
            if (inout_flag==0) then!读取入口断面
                ND=Com_element_main.Ca_cl(SEC_ID).CDS_in
                !-----------------------确定断面最低和最高点
                DO JJJ=1,ND!寻找高程点中的最大值
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
                    !若给定的水位比两个测点的高程均低，则循环变量J加1后继续执行循环
                    IF(ZW<=ZDI)THEN
                        CYCLE	 !CYCLE语句是终止本次循环，使循环变量加1再继续执行循环
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
            else !TU读取出口断面
                ND=Com_element_main.Ca_cl(SEC_ID).CDS_out
                !-----------------------确定断面最低和最高点
                DO JJJ=1,ND!寻找高程点中的最大值
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
                    !若给定的水位比两个测点的高程均低，则循环变量J加1后继续执行循环
                    IF(ZW<=ZDI)THEN
                        CYCLE	 !CYCLE语句是终止本次循环，使循环变量加1再继续执行循环
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
        !------------------------------------------#Tu修改到此为止
        !------------------------------------------
        CASE DEFAULT
    END SELECT
    !--------------------------------------------------
    RETURN
    END SUBROUTINE
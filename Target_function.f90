    !目标函数计算模块
    Subroutine Function_Target(GAV_Opendegree_F,ST_STATE_F,DS_STATE_F,Target1,Target2)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    REAL::Target1   !距离目标水位的最大偏差最小
    REAL::Target2   !水位平均偏差最小
    REAL::MAX_SWPC  !最大水位偏差
    INTEGER::i,j,Target1flag=0    !Target1flag:每小时水位变化是否超过允许值
    INTEGER::ADJUST_NUM !调控次数
    !REAL::Target2   !闸门调整次数

    REAL::ST_STATE_F((End_Time-Start_Time)/3600+1,2*Controlgate_NUM)!存储各节制闸每小时闸前水位、闸后水位
    REAL::DS_STATE_F((End_Time-Start_Time)/3600+1)!存储最下游边界水位

    REAL::GAV_Opendegree_F(1:Controlgate_NUM*Controltime)
    REAL::random_num

    !//////目标函数1计算过程///////平均值总是会出现问题，还是用最大偏差最小来处理
    Target1=0.0
    MAX_SWPC=0.0
    DO i=1,(End_Time-Start_Time)/3600+1
        DO j=1,Controlgate_NUM  !是否需要考虑第一个闸
            IF(ISNAN(ST_STATE_F(i,2*j-1)))THEN
                CALL Random_seed()  !系统根据日期和时间随机提供种子
                CALL RANDOM_NUMBER(random_num)
                Target1=Target1+1000+random_num
                EXIT
                !            ELSEIF(ST_STATE_F(i,2*j-1)>ST_CON_UP(j).OR.ST_STATE_F(i,2*j-1)<ST_CON_DOWN(J))THEN
                !                Target1=100000
                !                EXIT
            ELSE!执行这里的语句
                IF(ABS(ST_STATE_F(i,2*j-1)-ZQSW_TARGET(j))>=MAX_SWPC)THEN
                    MAX_SWPC=ABS(ST_STATE_F(i,2*j-1)-ZQSW_TARGET(j))!识别最大水位偏差然后输入MAX_SWPC
                ENDIF
                Target1=MAX_SWPC!然后把最大水位偏差输入Target1
            ENDIF
        ENDDO
    ENDDO
    DO i=1,(End_Time-Start_Time)/3600+1
        IF(ISNAN(DS_STATE_F(i)))THEN
            !Target1=100000
            CALL Random_seed()  !系统根据日期和时间随机提供种子
            CALL RANDOM_NUMBER(random_num)
            Target1=Target1+1000+random_num
            !write(*,*)"ISNAN(DS_STATE_F(i))"
            EXIT
            !        ELSEIF(DS_STATE_F(i)>DS_CON_UP.OR.DS_STATE_F(i)<DS_CON_DOWN)THEN
            !            Target1=100000
            !            EXIT
        ELSE!执行这里的
            IF(ABS(DS_STATE_F(i)-D_WL_TARGET)>=MAX_SWPC)THEN !D_WL_TARGET下游目标水位
                MAX_SWPC=ABS(DS_STATE_F(i)-D_WL_TARGET)
            ENDIF
            Target1=MAX_SWPC!事实上是节制闸闸前水位偏差以及最下游边界水位边界偏差的最大值
        ENDIF
    ENDDO

    !    IF(MAX_SWPC>=0.45.AND.MAX_SWPC<0.50)THEN
    !        Target1=10000.0
    !    ENDIF





    !   ////////////////////目标函数2计算过程 水位平均偏差最小
    Target2=0.0
    DO i=25,(End_Time-Start_Time)/3600+1
        DO j=1,Controlgate_NUM  !是否需要考虑第一个闸
            IF(ISNAN(ST_STATE_F(i,2*j-1)))THEN
                Target2=1000000
                EXIT
                !            ELSEIF(ST_STATE_F(i,2*j-1)>ST_CON_UP(j).OR.ST_STATE_F(i,2*j-1)<ST_CON_DOWN(J))THEN
                !                Target2=1000000
                !                EXIT
            ELSE
                Target2=Target2+ABS(ST_STATE_F(i,2*j-1)-ZQSW_TARGET(j))
            ENDIF
        ENDDO
    ENDDO
    DO i=25,(End_Time-Start_Time)/3600+1
        IF(ISNAN(DS_STATE_F(i)))THEN
            CALL Random_seed()  !系统根据日期和时间随机提供种子
            CALL RANDOM_NUMBER(random_num)
            Target2=1000000+random_num
            EXIT
            !        ELSEIF(DS_STATE_F(i)>DS_CON_UP.OR.DS_STATE_F(i)<DS_CON_DOWN)THEN
            !            Target2=1000000
            !            EXIT
        ELSE
            Target2=Target2+ABS(DS_STATE_F(i)-D_WL_TARGET)
        ENDIF
    ENDDO
    IF(ABS(Target2-1000000)>0.01)THEN
        !考虑第一个闸
        Target2=Target2/((End_Time-Start_Time)/3600+1-24)/(Controlgate_NUM+1)
        !不考虑第一个闸
        !Target2=Target2/((End_Time-Start_Time)/3600+1)/(Controlgate_NUM)
    ENDIF

    !闸门调控次数约束要求
    !    ADJUST_NUM=0
    !    DO i=1,156
    !        IF(ABS(GAV_Opendegree_F(i)>=0.03))THEN
    !            ADJUST_NUM=ADJUST_NUM+1
    !        ENDIF
    !    ENDDO
    !    IF(ADJUST_NUM>0.AND.ADJUST_NUM<26)THEN
    !        Target2=Target2+0
    !    ELSEIF(ADJUST_NUM>=26.AND.ADJUST_NUM<52)THEN
    !        Target2=Target2+0
    !    ELSEIF(ADJUST_NUM>=52.AND.ADJUST_NUM<78)THEN
    !        Target2=Target2+0
    !    ELSEIF(ADJUST_NUM>=78.AND.ADJUST_NUM<104)THEN
    !        Target2=Target2+0
    !    ELSEIF(ADJUST_NUM>=104.AND.ADJUST_NUM<130)THEN
    !        Target2=Target2+0
    !    ELSEIF(ADJUST_NUM>=130.AND.ADJUST_NUM<156)THEN
    !        Target2=Target2+5
    !    ENDIF

    !水位变幅约束要求

    
    
    
    IF(Target1/=100000.AND.Target2/=100000)THEN
        !每小时变幅不超过15cm
        DO j=1,Controlgate_NUM
            DO i=1,(End_Time-Start_Time)/3600
                IF((ST_STATE_F(i,2*j-1)-ST_STATE_F(i+1,2*j-1))>0.15)THEN
                    !Target1=100000
                    Target1=Target1+100*abs(ST_STATE_F(i,2*j-1)-ST_STATE_F(i+1,2*j-1))
                    Target1flag=1
                    !write(*,*)j,Start_Time
                    !write(*,*)ST_STATE_F(i,2*j-1)-ST_STATE_F(i+1,2*j-1)
                    !exit
                ENDIF
            ENDDO
        ENDDO
        DO i=1,(End_Time-Start_Time)/3600
            IF((DS_STATE_F(i)-DS_STATE_F(i+1))>0.15)THEN
                !Target1=100000
                Target1=Target1+100*abs(DS_STATE_F(i)-DS_STATE_F(i+1))
                exit
            ENDIF
        ENDDO
        !连续24小时变幅不超过30cm
        DO j=1,Controlgate_NUM
            DO i=1,(End_Time-Start_Time)/3600-24
                IF((ST_STATE_F(i,2*j-1)-ST_STATE_F(i+24,2*j-1))>0.3)THEN
                    !Target1=100000
                    Target1=Target1+100*abs(ST_STATE_F(i,2*j-1)-ST_STATE_F(i+24,2*j-1))
                    !Target1=Target1+1000
                    Target1flag=1
                    !exit
                ENDIF
            ENDDO
        ENDDO
        DO i=1,(End_Time-Start_Time)/3600-24
            IF((DS_STATE_F(i)-DS_STATE_F(i+24))>0.3)THEN
                !Target1=100000
                Target1=Target1+100*abs(DS_STATE_F(i)-DS_STATE_F(i+24))
                !Target1=Target1+1000
                !exit
            ENDIF
        ENDDO
    ENDIF
    IF(Target1<100.and.Target2<10) startflag=0
    IF(Target1<MIN_Target1) MIN_Target1=Target1
    IF(Target2<MIN_Target2) MIN_Target2=Target2
    !if(Target1>0.6)Target2=10000
    !write(*,*)Target1flag
    !write(*,*)Target1,Target2
    END
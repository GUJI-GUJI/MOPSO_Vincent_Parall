    !Ŀ�꺯������ģ��
    Subroutine Function_Target(GAV_Opendegree_F,ST_STATE_F,DS_STATE_F,Target1,Target2)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    REAL::Target1   !����Ŀ��ˮλ�����ƫ����С
    REAL::Target2   !ˮλƽ��ƫ����С
    REAL::MAX_SWPC  !���ˮλƫ��
    INTEGER::i,j,Target1flag=0    !Target1flag:ÿСʱˮλ�仯�Ƿ񳬹�����ֵ
    INTEGER::ADJUST_NUM !���ش���
    !REAL::Target2   !բ�ŵ�������

    REAL::ST_STATE_F((End_Time-Start_Time)/3600+1,2*Controlgate_NUM)!�洢������բÿСʱբǰˮλ��բ��ˮλ
    REAL::DS_STATE_F((End_Time-Start_Time)/3600+1)!�洢�����α߽�ˮλ

    REAL::GAV_Opendegree_F(1:Controlgate_NUM*Controltime)
    REAL::random_num

    !//////Ŀ�꺯��1�������///////ƽ��ֵ���ǻ�������⣬���������ƫ����С������
    Target1=0.0
    MAX_SWPC=0.0
    DO i=1,(End_Time-Start_Time)/3600+1
        DO j=1,Controlgate_NUM  !�Ƿ���Ҫ���ǵ�һ��բ
            IF(ISNAN(ST_STATE_F(i,2*j-1)))THEN
                CALL Random_seed()  !ϵͳ�������ں�ʱ������ṩ����
                CALL RANDOM_NUMBER(random_num)
                Target1=Target1+1000+random_num
                EXIT
                !            ELSEIF(ST_STATE_F(i,2*j-1)>ST_CON_UP(j).OR.ST_STATE_F(i,2*j-1)<ST_CON_DOWN(J))THEN
                !                Target1=100000
                !                EXIT
            ELSE!ִ����������
                IF(ABS(ST_STATE_F(i,2*j-1)-ZQSW_TARGET(j))>=MAX_SWPC)THEN
                    MAX_SWPC=ABS(ST_STATE_F(i,2*j-1)-ZQSW_TARGET(j))!ʶ�����ˮλƫ��Ȼ������MAX_SWPC
                ENDIF
                Target1=MAX_SWPC!Ȼ������ˮλƫ������Target1
            ENDIF
        ENDDO
    ENDDO
    DO i=1,(End_Time-Start_Time)/3600+1
        IF(ISNAN(DS_STATE_F(i)))THEN
            !Target1=100000
            CALL Random_seed()  !ϵͳ�������ں�ʱ������ṩ����
            CALL RANDOM_NUMBER(random_num)
            Target1=Target1+1000+random_num
            !write(*,*)"ISNAN(DS_STATE_F(i))"
            EXIT
            !        ELSEIF(DS_STATE_F(i)>DS_CON_UP.OR.DS_STATE_F(i)<DS_CON_DOWN)THEN
            !            Target1=100000
            !            EXIT
        ELSE!ִ�������
            IF(ABS(DS_STATE_F(i)-D_WL_TARGET)>=MAX_SWPC)THEN !D_WL_TARGET����Ŀ��ˮλ
                MAX_SWPC=ABS(DS_STATE_F(i)-D_WL_TARGET)
            ENDIF
            Target1=MAX_SWPC!��ʵ���ǽ���բբǰˮλƫ���Լ������α߽�ˮλ�߽�ƫ������ֵ
        ENDIF
    ENDDO

    !    IF(MAX_SWPC>=0.45.AND.MAX_SWPC<0.50)THEN
    !        Target1=10000.0
    !    ENDIF





    !   ////////////////////Ŀ�꺯��2������� ˮλƽ��ƫ����С
    Target2=0.0
    DO i=25,(End_Time-Start_Time)/3600+1
        DO j=1,Controlgate_NUM  !�Ƿ���Ҫ���ǵ�һ��բ
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
            CALL Random_seed()  !ϵͳ�������ں�ʱ������ṩ����
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
        !���ǵ�һ��բ
        Target2=Target2/((End_Time-Start_Time)/3600+1-24)/(Controlgate_NUM+1)
        !�����ǵ�һ��բ
        !Target2=Target2/((End_Time-Start_Time)/3600+1)/(Controlgate_NUM)
    ENDIF

    !բ�ŵ��ش���Լ��Ҫ��
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

    !ˮλ���Լ��Ҫ��

    
    
    
    IF(Target1/=100000.AND.Target2/=100000)THEN
        !ÿСʱ���������15cm
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
        !����24Сʱ���������30cm
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
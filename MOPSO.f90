    !	Copyright (c) 2016 Yousef Naranjani
    !	Email: yo.na1118@gmail.com
    !	Please read the Readme.txt file to get started
    PROGRAM MOPSOtools
    USE ProblemBank
    USE omp_lib
    USE GLOBAL_VARIABLE
    IMPLICIT NONE

    INTEGER:: i, j, gen, seed, time(8), funCnt = 0 ,kk    !endflag����5��Ŀ��仯С��һ��ֵ�����Ż�
    INTEGER:: fpLog, fpPF, fpPS, fpProgress, fpGenlog	! output file references
    REAL:: RAND ,xx, random_num  !random_num�����
    REAL:: dTarget1 ,dTarget2 ,tmpTarget1,tmpTarget2 = 0  !dTarget1,Target1ÿ���仯��   tmpTarget1,��ʱ�洢��һ�����
    real::mytemp
    real,allocatable::newopendegree(:)
    real,allocatable::iniinput(:)
    real::ini_lenth,ini_time
    integer::Cal_Time


    !ˮ�����½�����
    !TYPE(Boundary)::Boundary_in,Boundary_cal��ʱ���ȫ�ֱ���
    !TYPE(Config)::Config_main
    !TYPE(Com_element)::Com_element_main,Com_element_cal
    !TYPE(Initial_Result)::Initial_Result_main
    TYPE(Dynamic_Result)::Dynamic_Result_main
    !TYPE(Condition)::Condition_0

    INTEGER::k,gg
    REAL::A_n,B_n
    REAL(8)::start,finish



    integer :: maxArchive, gens, pops, NP, NF
    integer, allocatable :: N(:)
    real, allocatable :: lb(:), ub(:), lbFit(:), ubFit(:)
    logical :: generateReport, intermResults
    integer :: intermResultsInterval
    character(len = 30) :: probName
    integer :: indArchive, hyperLen
    integer*2, allocatable :: hyperSpace(:)
    real :: mutRate
    real, allocatable :: popArchive(:), fitArchive(:), pop(:), fit(:)
    integer, allocatable :: violArchive(:), violations(:), partPosArchive(:)
    real, allocatable :: popPbest(:), fitPbest(:)
    integer, allocatable :: violPbest(:)
    real, allocatable :: velocity(:)
    integer, allocatable :: lead(:)
    real :: t00, t01, t10, t11, t20, t21, t00_P, t01_P

    ! ����ӿں͹���ָ��
    abstract interface
    subroutine objective_function(fit, indv, viol, NF, NP)
    implicit none
    integer, intent(in) :: NF, NP
    real, intent(in) :: indv(NP)
    real, intent(out) :: fit(NF)
    integer, intent(out) :: viol
    end subroutine objective_function
    end interface

    procedure(objective_function), pointer :: f => null()
    ! f => MOP_CG_Vincent
    f => MOP_TU






    !INCLUDE 'in/MOP_CG_Vincent.INC'


    fpLog = 30
    fpPF = 31
    fpPS = 32
    fpProgress = 33
    fpGenlog = 34
    !open(unit=41,file='���ÿһ��բ�ſ���.csv')

    !CALL system('mkdir -p out/' // trim( probName ))

    !IF (generateReport) OPEN(UNIT = fpLog,FILE = ('out/log.txt'))
    !OPEN(UNIT = fpPF, FILE = ('out/PF.txt'))
    !OPEN(UNIT = fpPS, FILE = ('out/PS.txt'))
    !OPEN(UNIT = fpGenlog, FILE = ('out/genlog.txt'))
    !IF (intermResults) OPEN(UNIT = fpProgress,	FILE = ('out/progress.txt'))
    IF (generateReport) CALL DATE_AND_TIME(values=time)
    seed = 1000*time(7)+time(8)
    CALL SRAND(seed)

    !ˮ������ʼ����
    CALL Inputfile

    CALL Boundary_Meshing
    CALL Initial_Cal
    CALL Model_Meshing
    !ȫ�ֱ���������
    Boudarys_Number=Boundary_in%NUM
    Controlgate_NUMber=Com_element_main%Controlgate_NUM
    allocate(Opendegree_INIT(Boundary_in%Gate_NUM))
    Opendegree_INIT(:)=Boundary_in%Opendegree(:,1)

    allocate(iniinput(Boudarys_Number))
    iniinput(:)=0
    NP = Boudarys_Number
    Call Test_f(iniinput,ini_lenth,ini_time,NP)

    NP = (ini_time)/(Boundary_in.Step/60)

    IF(Boudarys_Number<NP+24)THEN
        write(*,*) "�����ı߽糤�ȹ��̣������ٸ���",(NP+24)*Boundary_in.Step,"��ı߽�����"
        stop
    ENDIF

    ! ��ʼ������
    maxArchive = 100
    gens = generations
    pops = populations
    !NP = 30
    NF = 2
    allocate(N(NF))
    N = (/30, 30/)
    allocate(lb(NP), ub(NP))
    lb = -0.3
    ub = 0.3
    generateReport = .FALSE.
    intermResults = .FALSE.
    intermResultsInterval = 20
    probName = "woxie MOP_CG_Vincent"
    allocate(lbFit(NF), ubFit(NF))
    lbFit = 0.0
    ubFit = 0.0
    indArchive = 0
    hyperLen = 1
    mutRate = 0.05
    allocate(popArchive(NP*maxArchive), fitArchive(NF*maxArchive), pop(NP*pops), fit(NF*pops))
    popArchive = 0.0
    fitArchive = 0.0
    pop = 0.0
    fit = 0.0
    allocate(violArchive(maxArchive), violations(pops), partPosArchive(maxArchive))
    violArchive = 0
    violations = 0
    partPosArchive = 0
    allocate(popPbest(NP*pops), fitPbest(NF*pops))
    popPbest = 0.0
    fitPbest = 0.0
    allocate(violPbest(pops))
    violPbest = 0
    allocate(velocity(NP*pops))
    velocity = 0.0
    allocate(lead(pops))



    ! ����hyperSpace�ĳ���
    do i = 1, NF
        hyperLen = hyperLen * N(i)
    end do
    allocate(hyperSpace(hyperLen))
    hyperSpace = 0


    !////////initialize and evaluate the population
    !CALL Random_seed()  !ϵͳ�������ں�ʱ������ṩ����!
    CALL init_random_seed()
    DO i = 1, pops
        DO j = 1, NP
            CALL RANDOM_NUMBER(random_num)!ÿ�θı������
            pop((i-1)*NP+j)=(random_num*(ub(j)-lb(j))+lb(j))
            !pop((i-1)*NP+j)=0
        END DO
    END DO
    !pop(15)=-3.33
    WRITE(*,*)"��ʼ��Ⱥ","������ϣ���"
    WRITE(*,*)

    !*********���м��㿪ʼʱ��********��
    CALL CPU_TIME(t00)
    !t00_P=omp_get_wtime()
    !Function evaluations


    CALL evaluate(fit, violations, pop, NF, NP, pops, funCnt, f)
    !mytemp=minVal(fit)
    !print*,mytemp
    !allocate(Ini_arriveTime(pops))
    !do i = 1,pops
    !    Ini_arriveTime(i)=fit(2*i)!������Ⱥ�м�������ĵ���ʱ�䣬��ȡ���е���Сֵ
    !end do
    !
    !Ini_CloseTime=Opendegree_INIT(1)/0.2
    !
    !write(*,*)
    !write(*,*)
    !write(*,*)"���������״����Ⱦ������ˮ�ڴ�Լ��Ҫ:",minval(Ini_arriveTime(:)),"����"
    !write(*,*)"Ӧ���رս���բ��Ҫ:",Ini_CloseTime,"����"
    !Cal_Time=int(minval(Ini_arriveTime(:))-Ini_CloseTime-20)
    !Cal_Time= Cal_Time*60/Boundary_in.step
    !Cal_Time=Cal_Time*Boundary_in.step!��λΪ��
    !write(*,*)"���ǵ���ȫ��������:",Cal_Time/60.0,"����  �ɹ����ɵ��ڽ���բ"
    !write(*,*)


    !Call Rebuild_by_Time(Cal_Time)!���ݼ���ʱ��ü�Inputfile
    !NP=Cal_Time/Boundary_in.step+1

    !	initializing the archive !��ʼ������
    CALL gBestUpdate(popArchive, fitArchive, violArchive, partPosArchive, indArchive, &
        pop, fit, violations, NP, NF, pops, maxArchive, &
        hyperSpace, hyperLen, lbFit, ubFit, N)
    !	initializing the local best ����ʼ����������ֵ
    popPbest(1:NP*pops) = pop(1:NP*pops)
    fitPbest(1:NF*pops) = fit(1:NF*pops)
    violPbest(1:pops) = violations(1:pops)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  THE FLIGHT CYCLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    gen = 0
    DO WHILE (endflag>0)
        IF(startflag) gen=0
        gen = gen + 1
        IF(gen==gens)then
            endflag=1
        endif

        WRITE(*,*)"--------��",gen,"��----------"
        !write(fpGenlog,*)
        !WRITE(fpGenlog,*)"--------��",gen,"��----------"
        !IF (generateReport) THEN
        !    WRITE(fpLog,"('Generation: ',I5,' began. Archive length: ', I5)") gen,indArchive
        !END IF

        ! choosing the leaders from Archive for velocity update
        CALL rouletteWheelSelection(lead, hyperSpace, hyperLen, indArchive, pops, partPosArchive)

        ! update the velocity here
        CALL velocityUpdate(velocity, pop, popPbest, popArchive, maxArchive, lead, pops, NP)

        ! Updating the population
        pop(1:NP*pops) = pop(1:NP*pops) + velocity(1:NP*pops)

        !����0-1��Ҫȡ��,������������բ��
        DO i = 1, pops
            DO j=1,NP
                IF(abs(pop((i-1)*NP+j))<0.03)THEN
                    pop((i-1)*NP+j)=0.0
                ENDIF
            ENDDO
        END DO

        ! Correcting the individuals that are out of bound
        CALL keepIn(pop, velocity, lb, ub, NP, pops)

        ! include mutation operator here
        CALL mutate(pop, pops, NP, gen, gens, mutRate, lb, ub)

        ! Correcting the individuals that are out of bound	 !������ٿ��Ƿ���ͻ�Ʊ߽��ֵ
        CALL keepIn(pop, velocity, lb, ub, NP, pops)

        !����0-1��Ҫȡ��,������������բ��
        DO i = 1, pops
            DO j=1,NP
                IF(abs(pop((i-1)*NP+j))<0.03)THEN
                    pop((i-1)*NP+j)=0.0
                ENDIF
            ENDDO
        END DO

        IF (generateReport) CALL CPU_TIME(t20)
        !	function evaluations
        CALL evaluate(fit, violations, pop, NF, NP, pops, funCnt, f)
        !mytemp=minVal(fit)
        !print*,'1minfit=',mytemp
        !IF (generateReport) THEN
        !    CALL CPU_TIME(t21)
        !    WRITE(fpLog,"('Function Eval: ',F10.2,' s')") t21-t20
        !END IF

        !	updating the archive

        CALL gBestUpdate(popArchive, fitArchive, violArchive, partPosArchive, indArchive, &
            pop, fit, violations, NP, NF, pops, maxArchive, &
            hyperSpace, hyperLen, lbFit, ubFit, N)

        !	updating the local best
        CALL pBestUpdate(popPbest, fitPbest, violPbest, pop, fit, violations, NP, NF, pops)

        !	write to progress
        IF ((intermResults.EQV. .TRUE.).AND.(MOD(gen,intermResultsInterval)==0)) THEN
            CALL writeTOfileProgress(fpProgress,	popArchive, fitArchive, &
                violArchive, gen, indArchive, NP, NF)
        END IF

        !IF (generateReport) THEN
        !    CALL CPU_TIME(t11)
        !    WRITE(fpLog,"('Generation: ',I5,' finished in ',F10.2,' s')") gen,t11-t10
        !END IF
        !CALL writeTOfile(fpGenlog, fitArchive, NF, indArchive)
        dTarget1 = abs(fit(1)-tmpTarget1)
        dTarget2 = abs(fit(2)-tmpTarget2)
        if(dTarget1<0.00002)then
            !if(dTarget2<0.01.and.fit(1)<10)then
            endflag=endflag-1
        else
            endflag=20
        endif

        IF(gen==gens)then
            endflag=0
        endif

        tmpTarget1=fit(1)
        tmpTarget2=fit(2)
        WRITE(*,*)'�����ɢ����Ϊ:',MIN_Target1,'��'
        WRITE(*,*)'��̵���ʱ��Ϊ:',MIN_Target2,'����'
        PRINT*,'endflag=',endflag
        WRITE(*,*)

        !WRITE(*,*)fit(1),MIN_Target1

    END DO
    !*********************************��
    !*********���м������ʱ��********��
    !*********************************!
    CALL CPU_TIME(t01)
    !t01_P =omp_get_wtime()

    !IF (generateReport) WRITE(fpLog,"('Total time: ',F10.2,' s')")t01-t00

    !CALL writeTOfile(fpPF, fitArchive, NF, indArchive)
    !CALL writeTOfile(fpPS, popArchive, NP, indArchive)

    !IF (generateReport) CLOSE(fpLog)
    !CLOSE(fpPF)
    !CLOSE(fpPS)
    !CLOSE(fpGenlog)
    !IF (intermResults) CLOSE(fpProgress)

    PRINT*,"===========����բ�Ż����̽�������ʼģ��Ӧ������ȫ����============"
    !PRINT*,'���м���ʱ��     ','Done! Time=',t01_P-t00_P,' ,',indArchive, ' solutions found'
    !PRINT*,'���м���ʱ��     ','Done! Time=',t01-t00,' ,',indArchive, ' solutions found'

    allocate(newopendegree(NP))
    Call Get_data(fitArchive,popArchive,NF,NP,indArchive,newopendegree)!���������������ǰ��


    Call Cal_After_Optimize(NP,newopendegree)!����ˮ��������

    !Call Cal_waterRetreating






    END PROGRAM MOPSOtools
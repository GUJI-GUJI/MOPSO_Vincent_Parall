    !	Copyright (c) 2016 Yousef Naranjani
    !	Email: yo.na1118@gmail.com
    !	Please read the Readme.txt file to get started
    PROGRAM MOPSOtools
    USE ProblemBank
    USE omp_lib
    USE GLOBAL_VARIABLE
    IMPLICIT NONE

    INTEGER:: i, j, gen, seed, time(8), funCnt = 0 ,kk    !endflag连续5代目标变化小于一定值结束优化
    INTEGER:: fpLog, fpPF, fpPS, fpProgress, fpGenlog	! output file references
    REAL:: RAND ,xx, random_num  !random_num随机数
    REAL:: dTarget1 ,dTarget2 ,tmpTarget1,tmpTarget2 = 0  !dTarget1,Target1每代变化量   tmpTarget1,临时存储上一代结果
    real::mytemp
    real,allocatable::newopendegree(:)
    real,allocatable::iniinput(:)
    real::ini_lenth,ini_time
    integer::Cal_Time


    !水动力新建变量
    !TYPE(Boundary)::Boundary_in,Boundary_cal暂时变成全局变量
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

    ! 抽象接口和过程指针
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
    !open(unit=41,file='输出每一次闸门开度.csv')

    !CALL system('mkdir -p out/' // trim( probName ))

    !IF (generateReport) OPEN(UNIT = fpLog,FILE = ('out/log.txt'))
    !OPEN(UNIT = fpPF, FILE = ('out/PF.txt'))
    !OPEN(UNIT = fpPS, FILE = ('out/PS.txt'))
    !OPEN(UNIT = fpGenlog, FILE = ('out/genlog.txt'))
    !IF (intermResults) OPEN(UNIT = fpProgress,	FILE = ('out/progress.txt'))
    IF (generateReport) CALL DATE_AND_TIME(values=time)
    seed = 1000*time(7)+time(8)
    CALL SRAND(seed)

    !水动力初始部分
    CALL Inputfile

    CALL Boundary_Meshing
    CALL Initial_Cal
    CALL Model_Meshing
    !全局变量先用着
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
        write(*,*) "给出的边界长度过短，请至少给出",(NP+24)*Boundary_in.Step,"秒的边界条件"
        stop
    ENDIF

    ! 初始化变量
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



    ! 计算hyperSpace的长度
    do i = 1, NF
        hyperLen = hyperLen * N(i)
    end do
    allocate(hyperSpace(hyperLen))
    hyperSpace = 0


    !////////initialize and evaluate the population
    !CALL Random_seed()  !系统根据日期和时间随机提供种子!
    CALL init_random_seed()
    DO i = 1, pops
        DO j = 1, NP
            CALL RANDOM_NUMBER(random_num)!每次改变随机数
            pop((i-1)*NP+j)=(random_num*(ub(j)-lb(j))+lb(j))
            !pop((i-1)*NP+j)=0
        END DO
    END DO
    !pop(15)=-3.33
    WRITE(*,*)"初始种群","生成完毕！！"
    WRITE(*,*)

    !*********并行计算开始时间********！
    CALL CPU_TIME(t00)
    !t00_P=omp_get_wtime()
    !Function evaluations


    CALL evaluate(fit, violations, pop, NF, NP, pops, funCnt, f)
    !mytemp=minVal(fit)
    !print*,mytemp
    !allocate(Ini_arriveTime(pops))
    !do i = 1,pops
    !    Ini_arriveTime(i)=fit(2*i)!整个种群中计算出来的到达时间，可取其中的最小值
    !end do
    !
    !Ini_CloseTime=Opendegree_INIT(1)/0.2
    !
    !write(*,*)
    !write(*,*)
    !write(*,*)"如果保持现状，污染到达退水口大约需要:",minval(Ini_arriveTime(:)),"分钟"
    !write(*,*)"应急关闭节制闸需要:",Ini_CloseTime,"分钟"
    !Cal_Time=int(minval(Ini_arriveTime(:))-Ini_CloseTime-20)
    !Cal_Time= Cal_Time*60/Boundary_in.step
    !Cal_Time=Cal_Time*Boundary_in.step!单位为秒
    !write(*,*)"考虑到安全余量，有:",Cal_Time/60.0,"分钟  可供自由调节节制闸"
    !write(*,*)


    !Call Rebuild_by_Time(Cal_Time)!根据计算时间裁剪Inputfile
    !NP=Cal_Time/Boundary_in.step+1

    !	initializing the archive !初始化档案
    CALL gBestUpdate(popArchive, fitArchive, violArchive, partPosArchive, indArchive, &
        pop, fit, violations, NP, NF, pops, maxArchive, &
        hyperSpace, hyperLen, lbFit, ubFit, N)
    !	initializing the local best ！初始化个体最优值
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

        WRITE(*,*)"--------第",gen,"代----------"
        !write(fpGenlog,*)
        !WRITE(fpGenlog,*)"--------第",gen,"代----------"
        !IF (generateReport) THEN
        !    WRITE(fpLog,"('Generation: ',I5,' began. Archive length: ', I5)") gen,indArchive
        !END IF

        ! choosing the leaders from Archive for velocity update
        CALL rouletteWheelSelection(lead, hyperSpace, hyperLen, indArchive, pops, partPosArchive)

        ! update the velocity here
        CALL velocityUpdate(velocity, pop, popPbest, popArchive, maxArchive, lead, pops, NP)

        ! Updating the population
        pop(1:NP*pops) = pop(1:NP*pops) + velocity(1:NP*pops)

        !对于0-1需要取整,并修正不调整闸门
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

        ! Correcting the individuals that are out of bound	 !交差后再看是否有突破边界的值
        CALL keepIn(pop, velocity, lb, ub, NP, pops)

        !对于0-1需要取整,并修正不调整闸门
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
        WRITE(*,*)'最短扩散距离为:',MIN_Target1,'米'
        WRITE(*,*)'最短到达时间为:',MIN_Target2,'分钟'
        PRINT*,'endflag=',endflag
        WRITE(*,*)

        !WRITE(*,*)fit(1),MIN_Target1

    END DO
    !*********************************！
    !*********并行计算结束时间********！
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

    PRINT*,"===========节制闸优化过程结束，开始模拟应急调控全过程============"
    !PRINT*,'并行计算时间     ','Done! Time=',t01_P-t00_P,' ,',indArchive, ' solutions found'
    !PRINT*,'串行计算时间     ','Done! Time=',t01-t00,' ,',indArchive, ' solutions found'

    allocate(newopendegree(NP))
    Call Get_data(fitArchive,popArchive,NF,NP,indArchive,newopendegree)!获得排序后的帕累托前沿


    Call Cal_After_Optimize(NP,newopendegree)!进行水动力计算

    !Call Cal_waterRetreating






    END PROGRAM MOPSOtools
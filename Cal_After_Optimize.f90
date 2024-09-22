    subroutine Cal_After_Optimize(NP,oldOpendegree)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    INTEGER, INTENT(IN)::NP
    Real::oldOpendegree(NP)
    Integer::add_CalTime
    Integer::Endtime!关闸后的时间
    Integer::i,j
    real::Up_tmp,Down_tmp
    Integer::changetime

    add_CalTime=int(oldOpendegree(size(oldOpendegree))/0.2*60/boundary_in.step)+1!应急关闸所需的时间
    add_CalTime=add_CalTime+12

    Boundary_2in=Boundary_in
    Boundary_2in.End_Time=Boundary_2in.Start_Time+(NP+add_CalTime-1)*Boundary_2in.Step
    CALL Boundary_Auto_Update(Boundary_2in,"NUM")
    CALL Boundary_Auto_Update(Boundary_2in,"Gate_NUM")
    CALL Boundary_Auto_Update(Boundary_2in,"SideOutFlow_NUM")


    !获取Opendegree
    DO i=1,Com_element_main.Controlgate_NUM
        Boundary_2in.Opendegree(i,1:NP/Com_element_main.Controlgate_NUM)=oldOpendegree(1+(i-1)*NP/Com_element_main.Controlgate_NUM:i*NP/Com_element_main.Controlgate_NUM)
    ENDDO

    !后面紧急关闸
    Boundary_2in.Opendegree(Com_element_main.Controlgate_NUM,NP+1:size(Boundary_2in.Opendegree,2))=0!注：紧急关闸过后的水位可能暴涨
    !do i=1,Com_element_main.Controlgate_NUM-1
    !    Boundary_2in.Opendegree(i,NP/Com_element_main.Controlgate_NUM+1:size(Boundary_2in.Opendegree,2))=Boundary_2in.Opendegree(i,NP/Com_element_main.Controlgate_NUM)
    !enddo
    do i = 1,add_CalTime-1-12
        Boundary_2in.Opendegree(Com_element_main.Controlgate_NUM,NP+i)=oldOpendegree(size(oldOpendegree))-1*i!每五分钟关1m
    end do


    DO i=1,NP+add_CalTime !读边界条件
        Boundary_2in.DownZ(i)=Boundary_in.downZ(1)
    END DO

    DO i=1,NP !读边界条件
        Boundary_2in.UpQ(i)=Boundary_in.UpQ(1)
    END DO
    do i = Np+1,NP+add_CalTime
        Boundary_2in.UpQ(i)= MAX(Boundary_2in.UpQ(i-1) - inGate_changelimit*Boundary_in.Step/60, 0.0)
    end do


    Boundary_2in.OutDischarge(:,:)=0!!!!!!!!!!!!
    !Do i = 1,add_CalTime
    !    !Boundary_2in.OutDischarge(2,NP+i)=10*i!每五分钟的流量变化
    !    Boundary_2in.OutDischarge(:,NP/Com_element_main.Controlgate_NUM+i)=0!每五分钟的流量变化,暂时先这样
    !end do



    changetime = int((outflow_maxQ/(outflow_changelimit*(Boundary_in.Step/60))))
    !开退水口
    do i = 1,changetime
        Boundary_2in.OutDischarge(outflow_index,NP+i-3)=outflow_changelimit*Boundary_in.Step/60*i
    end do


    Boundary_2in.OutDischarge(outflow_index,NP+changetime-2:)=outflow_maxQ



    Call Boundary_Meshing_fun(Boundary_2in,Boundary_2cal)

    Endtime = NP+add_CalTime-12
    Endtime = Boundary_2in.Step/Boundary_2cal.Step*(Endtime-1)+1

    CALL Dynamic_Cal2(Boundary_2cal,Up_tmp,down_tmp,Endtime)!后面两个凑数的


    !write(*,*)"退水闸调控过程从",Boundary_2in.NUM,"个边界开始计算"
    !write(*,*)"已经关闸了",add_CalTime,"个边界时刻"





    end subroutine
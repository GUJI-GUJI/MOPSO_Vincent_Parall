    subroutine Cal_After_Optimize(NP,oldOpendegree)
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    INTEGER, INTENT(IN)::NP
    Real::oldOpendegree(NP)
    Integer::add_CalTime
    Integer::Endtime!��բ���ʱ��
    Integer::i,j
    real::Up_tmp,Down_tmp
    Integer::changetime

    add_CalTime=int(oldOpendegree(size(oldOpendegree))/0.2*60/boundary_in.step)+1!Ӧ����բ�����ʱ��
    add_CalTime=add_CalTime+12

    Boundary_2in=Boundary_in
    Boundary_2in.End_Time=Boundary_2in.Start_Time+(NP+add_CalTime-1)*Boundary_2in.Step
    CALL Boundary_Auto_Update(Boundary_2in,"NUM")
    CALL Boundary_Auto_Update(Boundary_2in,"Gate_NUM")
    CALL Boundary_Auto_Update(Boundary_2in,"SideOutFlow_NUM")


    !��ȡOpendegree
    DO i=1,Com_element_main.Controlgate_NUM
        Boundary_2in.Opendegree(i,1:NP/Com_element_main.Controlgate_NUM)=oldOpendegree(1+(i-1)*NP/Com_element_main.Controlgate_NUM:i*NP/Com_element_main.Controlgate_NUM)
    ENDDO

    !���������բ
    Boundary_2in.Opendegree(Com_element_main.Controlgate_NUM,NP+1:size(Boundary_2in.Opendegree,2))=0!ע��������բ�����ˮλ���ܱ���
    !do i=1,Com_element_main.Controlgate_NUM-1
    !    Boundary_2in.Opendegree(i,NP/Com_element_main.Controlgate_NUM+1:size(Boundary_2in.Opendegree,2))=Boundary_2in.Opendegree(i,NP/Com_element_main.Controlgate_NUM)
    !enddo
    do i = 1,add_CalTime-1-12
        Boundary_2in.Opendegree(Com_element_main.Controlgate_NUM,NP+i)=oldOpendegree(size(oldOpendegree))-1*i!ÿ����ӹ�1m
    end do


    DO i=1,NP+add_CalTime !���߽�����
        Boundary_2in.DownZ(i)=Boundary_in.downZ(1)
    END DO

    DO i=1,NP !���߽�����
        Boundary_2in.UpQ(i)=Boundary_in.UpQ(1)
    END DO
    do i = Np+1,NP+add_CalTime
        Boundary_2in.UpQ(i)= MAX(Boundary_2in.UpQ(i-1) - inGate_changelimit*Boundary_in.Step/60, 0.0)
    end do


    Boundary_2in.OutDischarge(:,:)=0!!!!!!!!!!!!
    !Do i = 1,add_CalTime
    !    !Boundary_2in.OutDischarge(2,NP+i)=10*i!ÿ����ӵ������仯
    !    Boundary_2in.OutDischarge(:,NP/Com_element_main.Controlgate_NUM+i)=0!ÿ����ӵ������仯,��ʱ������
    !end do



    changetime = int((outflow_maxQ/(outflow_changelimit*(Boundary_in.Step/60))))
    !����ˮ��
    do i = 1,changetime
        Boundary_2in.OutDischarge(outflow_index,NP+i-3)=outflow_changelimit*Boundary_in.Step/60*i
    end do


    Boundary_2in.OutDischarge(outflow_index,NP+changetime-2:)=outflow_maxQ



    Call Boundary_Meshing_fun(Boundary_2in,Boundary_2cal)

    Endtime = NP+add_CalTime-12
    Endtime = Boundary_2in.Step/Boundary_2cal.Step*(Endtime-1)+1

    CALL Dynamic_Cal2(Boundary_2cal,Up_tmp,down_tmp,Endtime)!��������������


    !write(*,*)"��ˮբ���ع��̴�",Boundary_2in.NUM,"���߽翪ʼ����"
    !write(*,*)"�Ѿ���բ��",add_CalTime,"���߽�ʱ��"





    end subroutine
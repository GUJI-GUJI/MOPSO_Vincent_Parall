    MODULE GLOBAL_VARIABLE  !�������
    IMPLICIT NONE
    !����ȫ�ֱ���

    !V2.0V2.0V2.0
    !V2.0V2.0V2.0
    !V2.0V2.0V2.0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!��������!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    REAL,PARAMETER :: PI=3.14159265

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!���̲�������!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !////��������//////////
    Type::Canal
        INTEGER::Serial_NUM                 !�ֲ����
        INTEGER::irr_Serial_NUM=0           !��ʾ�������ǲ����������еĵڼ���
        CHARACTER(len=12)::Name             !�ڵ�����
        REAL::Incoord                       !������
        REAL::Outcoord                      !�������
        REAL::BottomWidth                   !�׿�
        REAL::SideSlopes                    !����ϵ��
        REAL::InInvertElev                  !��ڵ׸߳�
        REAL::OutInvertElev                 !���ڵ׸߳�
        REAL::InTopElev                     !��ڶ��߳�
        REAL::OutTopElev                    !���ڶ��߳�
        REAL::Roughness                     !��������
        INTEGER :: CroSecTYPE               !�������ͣ�0 ���Σ����Σ���1Բ�Σ�2���������
        INTEGER  :: CDS_in                  !#TU��ڶ�������
        INTEGER  :: CDS_out                 !#TU��ڶ�������
        REAL,DIMENSION(2, 3000) :: GCD_in   !#TU��ڸ̶߳�
        REAL,DIMENSION(2, 3000) :: GCD_out  !#TU���ڸ̶߳�
    End type Canal

    !/////���������////////////////
    Type::Areachange
        INTEGER::Serial_NUM         !�ֲ����
        CHARACTER(len=12)::Name     !�ڵ�����
        REAL::Incoord               !������
        REAL::Outcoord              !�������
        REAL::Locallosscoeff        !�ֲ�ˮ����ʧϵ��
        REAL::Additionlosscoeff     !������ʧϵ��
        REAL::SideSlopes            !����ϵ��
    End type Areachange

    !/////����������////////////////
    Type::Invertedsiphon
        INTEGER::Serial_NUM     !�ֲ����
        CHARACTER(len=12)::Name !�ڵ�����
        REAL::Incoord           !������
        REAL::Outcoord          !�������
        INTEGER::ParallelNum    !���Ÿ���
        REAL::BottomWidth       !���
        REAL::Height            !�߶�
        REAL::Length            !����
        REAL::InInvertElev      !��ڵ׸߳�
        REAL::OutInvertElev     !���ڵ׸߳�
        REAL::Roughness         !����������
        REAL::Insi_coeff        !�ֲ�ˮ����ʧϵ��
    End type Invertedsiphon

    !////�Ŷ�����//////////
    Type::Bridge
        INTEGER::Serial_NUM     !�ֲ����
        CHARACTER(len=12)::Name !�ڵ�����
        REAL::Incoord           !������
        REAL::Outcoord          !�������
        REAL::BottomWidth       !�׿�
        REAL::SideSlopes        !����ϵ��
        REAL::InInvertElev      !��ڵ׸߳�
        REAL::OutInvertElev     !���ڵ׸߳�
        REAL::Roughness         !��������
        REAL::HinderWidth       !����ˮ���
    End type Bridge

    !/////����բ����////////////////
    Type::Controlgate
        INTEGER::Serial_NUM     !�ֲ����
        CHARACTER(len=12)::Name !�ڵ�����
        REAL::Incoord           !������
        REAL::Outcoord          !�������
        REAL::Initial_gate_Zup  !��ʼˮλ
        INTEGER::ParallelNum    !���Ÿ���
        REAL::SingleWidth       !���
        REAL::InvertElev        !�׸߳�
        REAL::Opendegree        !բ�ſ���
        INTEGER:: RUNNum        !��������
        REAL::G_A,G_B,G_C       !��բ����ϵ����բ�ſ���֮��Ĺ�ϵ����
        REAL::Discoeff          !բ�Ź���ϵ��
        REAL::Lenth
    End type Controlgate

    !/////��ˮ������////////////////
    Type::SideOutFlow
        INTEGER::Serial_NUM     !�ֲ����
        CHARACTER(len=12)::Name !�ڵ�����
        REAL::Incoord           !������
        REAL::Outcoord          !�������
        REAL::OutDischarge      !��ˮ����
    End type SideOutFlow

    !/////��վ����////////////////
    Type::Pumpingstation
        INTEGER::Serial_NUM                 !�ֲ����
        CHARACTER(len=12)::Name             !�ڵ�����
        REAL::Incoord                       !������
        REAL::Outcoord                      !�������
        REAL::Initial_pump_Zup              !��ʼˮλ
        REAL::InvertElev                    !�׸߳�
        INTEGER::UNIT_NUM                   !��վ�������
        INTEGER,ALLOCATABLE::Pump_Type(:)   !��վ��ˮ�����͵ı��
        REAL,ALLOCATABLE::Coefficient_A(:)  !������ϵ������
        REAL,ALLOCATABLE::Coefficient_B(:)  !һ����ϵ������
        REAL,ALLOCATABLE::Coefficient_C(:)  !����������
    End type Pumpingstation

    !/////ˮ������////////////////
    Type::Pump
        INTEGER::Serial_NUM                     !���ͱ��
        CHARACTER(len=12)::Name                 !ˮ������
        INTEGER::Coefficient_NUM                !ת�Ǹ���
        REAL,ALLOCATABLE::Bladeangle_Serial(:)  !ת������
        REAL,ALLOCATABLE::Coefficient_A(:)      !������ϵ������
        REAL,ALLOCATABLE::Coefficient_B(:)      !һ����ϵ������
        REAL,ALLOCATABLE::Coefficient_C(:)      !����������
    END TYPE Pump

    !/////ת�ǹ���ʱ������//////////
    Type::Bladeangle_TimeSerial
        REAL,ALLOCATABLE::Each_PS(:,:)   !��վ��ÿ���õ�ת�ǹ���   ��ת�ǣ�ʱ�䣩
    END TYPE Bladeangle_TimeSerial

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!��������!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !�߽�����
    TYPE::Boundary
        CHARACTER(LEN=50)::name             !����
        CHARACTER(LEN=50)::Filename         !�ļ���ַ
        INTEGER::Start_Time                 !��ʼʱ��
        INTEGER::End_Time                   !����ʱ�䣨�����ʾ��
        INTEGER::Step                       !�߽���������
        INTEGER::Upsign                     !�߽��������Ʒ���0-������1-ˮλ��2-ˮλ~������ϵ
        INTEGER::Downsign                   !�߽��������Ʒ���0-������1-ˮλ��2-ˮλ~������ϵ
        REAL,ALLOCATABLE::UpQ(:)            !������������
        REAL,ALLOCATABLE::UpZ(:)            !����ˮλ����
        REAL,ALLOCATABLE::DownQ(:)          !������������
        REAL,ALLOCATABLE::DownZ(:)          !����ˮλ����
        REAL,ALLOCATABLE::Opendegree(:,:)   !բ�ſ��ȹ���
        REAL,ALLOCATABLE::OutDischarge(:,:) !��ˮ�ڳ�������
        Type(Bladeangle_TimeSerial),&
            ALLOCATABLE::Bladeangle_cl(:)   !�������Ǳ�վ���

        LOGICAL::Auto_Update                !�Ƿ��Զ��������в��������Զ���������check�н��и��£�
        INTEGER::NUM                        !�߽����и������߽������еĹ۲����ݣ�
        INTEGER::Gate_NUM                   !բ������
        INTEGER::PumpST_NUM                !��վ����
        INTEGER::SideOutFlow_NUM            !��ˮ������
    END TYPE Boundary

    !������Ϣ
    TYPE::Com_element
        CHARACTER(LEN=50)::name                         !����
        CHARACTER(LEN=50)::Filename                     !�ļ���ַ
        INTEGER::Element_NUM                            !����Ԫ���ĸ���
        INTEGER,ALLOCATABLE::Ele_Rel(:,:)               !element_relation�洢ϵͳ��š����ʹ��롢�ֲ����
        INTEGER::Canal_NUM                              !��������
        Type(Canal),ALLOCATABLE::Ca_cl(:)               !canal class�趨��������
        INTEGER::Areachange_NUM                         !����θ���
        Type(Areachange),ALLOCATABLE::Arch_cl(:)        !areachange class �趨���������
        INTEGER::Invertedsiphon_NUM                     !����������
        Type(Invertedsiphon),ALLOCATABLE::Insi_cl(:)    !Invertedsiphon class �趨����������
        INTEGER::Controlgate_NUM                        !բ�Ÿ���
        Type(Controlgate),ALLOCATABLE::Cogate_cl(:)     !Controlgate class �趨����բ����
        INTEGER::SideOutFlow_NUM                        !��ˮ�ڸ���
        Type(SideOutFlow),ALLOCATABLE::SideOF_cl(:)     !SideOutFlow class �趨��ˮ������
        INTEGER::Bridge_NUM                             !�Ŷո���
        Type(Bridge),ALLOCATABLE::Brid_cl(:)            !Bridge class�趨��������
        INTEGER::PumpST_NUM                             !��վ����
        Type(Pumpingstation),ALLOCATABLE::PumpST_cl(:)  !Pumpingstation class �趨��վ����
        INTEGER::Pump_type_NUM                          !ˮ�����͸���
        Type(Pump),ALLOCATABLE::Pump_cl(:)              !Pumpingstation class �趨ˮ������

        LOGICAL::Auto_Update                            !�Ƿ��Զ��������в��������Զ���������check�н��и��£�
        INTEGER,ALLOCATABLE::n(:)                       !ÿһ��Ԫ���ķֶ����
        INTEGER,ALLOCATABLE::n_canal(:)                 !ÿ�������ķֶ����
        INTEGER,ALLOCATABLE::SG_POSTION(:)              !����բλ��
        INTEGER,ALLOCATABLE::DM_POSTION(:)              !����բ+��β����λ��
        INTEGER,ALLOCATABLE::ELE_POSTION(:)             !�������Ŷ������Ԫ��λ��
    END TYPE Com_element

    !��Ӧʱ�̵ĸ�����ֵ
    TYPE::Condition
        CHARACTER(LEN=50)::Comname          !��Ӧ������Ϣ��
        CHARACTER(LEN=50)::Boundaryname     !��Ӧ�߽�������Ϣ��
        REAL::time                          !ʱ��
        REAL::h(10000)                      !��ʼˮ�ԭwaterdepth
        REAL::Q(10000)                      !��ʼ������ԭwaterdepth
        REAL,ALLOCATABLE::ZZ(:,:)           !ˮλ
        REAL,ALLOCATABLE::QQ(:,:)           !����
        REAL,ALLOCATABLE::ELE_Volume(:,:)   !����
        REAL,ALLOCATABLE::SS(:)             !׷��ϵ��
        REAL,ALLOCATABLE::TT(:)             !׷��ϵ��
        REAL,ALLOCATABLE::PP(:)             !׷��ϵ��
        REAL,ALLOCATABLE::VV(:)             !׷��ϵ��

        LOGICAL::Auto_Update                !�Ƿ��Զ��������в���
        REAL,ALLOCATABLE::Pool_Volume(:,:)  !��������
        REAL,ALLOCATABLE::Total_Volume(:)   !������
    END TYPE Condition

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!ģ�Ͳ�������!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !�����������������ģ���ȶ��Բ����Ͳ��ʡ��ֲ�ˮ����ʧ��
    TYPE::Coeff
        REAL::Theta                         !�ռ�Ȩ��ϵ��(��������)
        REAL::Phi                           !ʱ��Ȩ��ϵ��
        REAL::Psi                           !�ռ�Ȩ��ϵ��(��������)
        REAL::PsiR                          !�ռ�Ȩ��ϵ��(Ħ���¶�)
        REAL::Beta                          !���涯��У��ϵ��
        REAL::GRAV                          !�������ٶ�
        REAL::RoughnessCoeff                !���� ��������������ͬһ������ϵ����
        REAL,ALLOCATABLE::Insi_coeff(:,:)   !��������ʧϵ����Ŀǰ����2����ʧϵ����������Ҫ���е�����
        REAL::Discoeff_GATE                 !բ�Ź���ϵ��
    END TYPE Coeff

    !ͬ������
    TYPE::Assimilate
        INTEGER::NUMMember      !�����г�Ա����
        INTEGER::Output_sign    !���������Ʒ���0-ֻ�������բ�������1-�����Ԫ��������
        REAL::m_coeff           !ĩ��բ�Ź�բ����ϵ������
        INTEGER::Control_sign   !ͬ��ģ�ͺ�ģ��ģ�͵Ŀ��Ʒ���ģ��Ϊ0��ͬ��Ϊ1
        INTEGER::Statetype      !״̬�����ͣ�0-���ʣ�1-����2-����+����
        INTEGER::Obsertype      !�۲������ͣ�0-ˮλ��1-������2-ˮλ+������
        REAL::PPQ               !ˮλ~������ϵ����ϵ��
        REAL::VVQ               !ˮλ~������ϵ����ϵ��
    END TYPE Assimilate

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!ģ�Ͳ���!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !ģ��������
    TYPE::Config
        CHARACTER(LEN=50)::Filename         !�ļ���ַ
        TYPE(Coeff)::Coeff_set
        TYPE(Assimilate)::Assimilate_set
    END TYPE Config

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!������!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !��ʼ�㶨��������
    TYPE::Initial_Result
        CHARACTER(LEN=50)::Filename         !�ļ���ַ
        CHARACTER(LEN=50)::Comname          !��Ӧ������Ϣ��
        CHARACTER(LEN=50)::Boundaryname     !��Ӧ�߽�������Ϣ��
        REAL,ALLOCATABLE::h(:)              !ˮ��
        REAL,ALLOCATABLE::Q(:)              !����
        REAL,ALLOCATABLE::Z(:)              !ˮλ
        REAL,ALLOCATABLE::vc(:)             !����
        REAL,ALLOCATABLE::L(:)              !���
        REAL,ALLOCATABLE::Vol(:)            !����
        REAL,ALLOCATABLE::Pool_Volume(:)    !��������

        LOGICAL::Auto_Update                !�Ƿ��Զ��������в���
        INTEGER::NUM                        !Ԫ������
    END TYPE Initial_Result

    !�Ǻ㶨�����
    TYPE::Dynamic_Result
        CHARACTER(LEN=50)::Filename         !�ļ���ַ
        CHARACTER(LEN=50)::Comname          !��Ӧ������Ϣ��
        CHARACTER(LEN=50)::Boundaryname     !��Ӧ�߽�������Ϣ��
        INTEGER::DeltaT                     !ʱ�䲽��
        REAL,ALLOCATABLE::Zup(:,:)          !բǰˮλ
        REAL,ALLOCATABLE::Zdown(:,:)        !բ��ˮλ
        REAL,ALLOCATABLE::Q(:,:)            !��բ����
        REAL,ALLOCATABLE::Zbd(:,:)          !ĩ��ˮλ
        REAL,ALLOCATABLE::Qbd(:,:)          !ĩ������
        REAL,ALLOCATABLE::ZZ(:,:)           !������ˮλ
        REAL,ALLOCATABLE::QQ(:,:)           !����������
        REAL,ALLOCATABLE::Pool_Volume(:,:)    !��������
        REAL,ALLOCATABLE::Total_Volume(:)   !������

        INTEGER::NUM
    END TYPE Dynamic_Result

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    INTEGER::Out_DeltaT !����ʱ�䲽����������������ʱ�䲽��
    INTEGER::Output_sign    !���������Ʒ���0-ֻ�������բ�������1-�����Ԫ��������

    INTEGER::PumpST_NUM    !��������������θ�����������������բ�Ÿ���,��ˮ�ڸ���,�Ŷո�������վ����
    INTEGER::Pump_type_NUM !ˮ�����͸���
    !//////////////////////////////////////////////
    Type(Pumpingstation),ALLOCATABLE::PumpST_cl(:)    !Pumpingstation class �趨��վ����
    Type(Pumpingstation),ALLOCATABLE::Ini_PumpST_cl(:)
    !//////////////////////////////////////////////
    Type(Pump),ALLOCATABLE::Pump_cl(:)    !Pumpingstation class �趨ˮ������

    REAL,ALLOCATABLE::C_Bladeangle(:,:,:)   !�������ʹ��բ�ſ���
    Type(Bladeangle_TimeSerial),ALLOCATABLE::Bladeangle_cl(:)       !�������Ǳ�վ���
    Type(Bladeangle_TimeSerial),ALLOCATABLE::C_Bladeangle_cl(:)       !�������Ǳ�վ���






    !New����MPC
    REAL::MIN_Target1=10000,MIN_Target2=10000
    REAL,ALLOCATABLE::ZQSW_TARGET(:)    !����բբǰĿ��ˮλ���м�������բ���м���
    REAL::D_WL_TARGET   !�����α߽紦��Ŀ��ˮλ
    REAL,ALLOCATABLE::ST_CON_UP(:)  !ˮλԼ������
    REAL,ALLOCATABLE::ST_CON_DOWN(:)    !ˮλԼ������
    REAL::DS_CON_UP,DS_CON_DOWN !����ˮλ�߽磬���ޡ�����
    REAL,ALLOCATABLE::h_new_INIT(:),QQ_INIT(:,:)    !��ʼˮλ����ʼ����
    !REAL,ALLOCATABLE::Opendegree_INIT(:,:)  !��ʼ����


    REAL,ALLOCATABLE::BestOpendegree(:,:)
    REAL,ALLOCATABLE::BestST_STATE(:,:)
    REAL,ALLOCATABLE::BestQControlGate(:,:)
    REAL,ALLOCATABLE::BestVolume(:,:)




    !//////���ѭ������
    INTEGER::Cycle_NUM  !ÿһ�����������ѭ����
    !///////���Է������е�δ֪��

    !/////��������
    INTEGER::endflag=20  !5�������������������

    !/////��ʼ����
    INTEGER::startflag=1

    !/////����������
    INTEGER::bestoutputflag=1





    !temp�������ú��񲻺�
    integer::Boudarys_Number
    integer::Controlgate_NUMber
    real,allocatable::Opendegree_INIT(:)!��ȡ���ȵĵ�һ��
    TYPE(Boundary)::Boundary_in,Boundary_cal,Boundary_2in,Boundary_2cal
    TYPE(Config)::Config_main
    TYPE(Com_element)::Com_element_main,Com_element_cal
    TYPE(Initial_Result)::Initial_Result_main
    TYPE(Condition)::Condition_0
    
    Integer::Function_index=0!evaluate���������д���
    Real,allocatable::Ini_ArriveTime(:)!Ԥ�����Ⱦ����ʱ��
    Real::Ini_CloseTime!Ԥ���բ�Źر��������ʱ��
    Integer::total_BoundaryNum=97
    
    Integer::firstgate_Index
    Integer::firstGateIndex
    Integer::outflow_index
    Integer::arrivetime_weight
    Real::pollut_location!��ȾԴ�������ξ���
    Real::pollut_weight
    Real::pollut_ConveyDistance
    Real::outflow_maxQ
    Real::outflow_changelimit
    Real::inGate_changelimit
    INTEGER::generations
    INTEGER::populations
    !Integer::maxArchive,gens,pops,NP,NF
    End


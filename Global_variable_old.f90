MODULE GLOBAL_VARIABLE  !�������
    IMPLICIT NONE
    !����ȫ�ֱ���
    !�����ļ�
    INTEGER::Start_Time,End_Time   !��ʼʱ�䣬����ʱ�䣨�����ʾ��
    INTEGER::Boundary_Step,Boundary_NUM !�߽������������߽����и������߽������еĹ۲����ݣ�
    INTEGER::Controltime    !���ش���
    INTEGER::DeltaT,Calculate_NUM,Out_DeltaT !����ʱ�䲽����������������ʱ�䲽��
    REAL::DeltaX    !�ռ䲽��
    REAL,ALLOCATABLE::delX(:) !�洢ÿһ�������ķֶ���
    INTEGER::c_element    !��������ռ䲽����������Ԫ���ĸ���
    INTEGER,ALLOCATABLE::nnn(:)   !ÿһ��Ԫ���ķֶ����
    INTEGER,ALLOCATABLE::n_canal(:) !ÿ�������ķֶ����
    REAL,ALLOCATABLE::R_Upstream(:),R_Downstream(:) !���α߽����������α߽�����
    REAL,ALLOCATABLE::Upstream(:),Downstream(:) !���α߽����������α߽�����
    INTEGER::Element_NUM    !����Ԫ���ĸ���
    INTEGER::Canal_NUM,Areachange_NUM,Invertedsiphon_NUM,Controlgate_NUM,SideOutFlow_NUM,Bridge_NUM,PumpST_NUM    !��������������θ�����������������բ�Ÿ���,��ˮ�ڸ���,�Ŷո�������վ����
    INTEGER::Pump_type_NUM !ˮ�����͸���
    INTEGER::Canal_New_NUM  !ʵ�ʲ���������εĸ���
    !��ϵͳ��ű�ʶ����Ԫ�������������е�˳�򡣾ֲ����Ϊ��ͬ��Ԫ���е���š����ʹ����ʾ����Ԫ��������
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation�洢ϵͳ��š����ʹ��롢�ֲ����
    INTEGER,ALLOCATABLE::Ele_Rel_New(:,:)   !element_relation_new�洢����ʱ�䲽���Ϳռ䲽��֮���ϵͳ��š����ʹ��롢�ֲ����
    REAL:: Initial_Discharge    !���α߽�Ҳѡ��ˮλʱ�ĳ�ʼ����
    REAL:: Initial_gate_Zup     !��ʼբǰˮλ
    REAL:: Initial_pump_Zup     !��ʼվǰˮλ
    REAL:: D_elevation  !���α߽�׸߳�
    REAL::Waterdepth(10000,2),Waterdepth_flag
    REAL,ALLOCATABLE::Dynamic_Result_Zup(:,:),Dynamic_Result_Zdown(:,:),Dynamic_Result_Q(:,:),Dynamic_Result_Zbd(:,:),Dynamic_Result_Qbd(:,:),Dynamic_Result_ZZ(:,:),Dynamic_Result_QQ(:,:)
    REAL::MIN_Target1=10000,MIN_Target2=10000
    
    !/////�����Ϣ
    REAL,ALLOCATABLE::h(:),Q(:),Z(:),vc(:),L(:),Vol(:)
    
    !////��������//////////
    Type::Canal
        INTEGER::Serial_NUM !�ֲ����
        CHARACTER(len=12)::Name !�ڵ�����
        REAL::Incoord !������
        REAL::Outcoord !�������
        REAL::BottomWidth   !�׿�
        REAL::SideSlopes    !����ϵ��
        REAL::InInvertElev  !��ڵ׸߳�
        REAL::OutInvertElev  !���ڵ׸߳�
        REAL::InTopElev !��ڶ��߳�
        REAL::OutTopElev    !���ڶ��߳�
        REAL::Roughness !��������
        INTEGER :: CroSecTYPE  !�������ͣ�0 ���Σ����Σ���1Բ�Σ�2��������� !wjb0414
        INTEGER  :: CDS  !��������, ����������� ID_IR_CroSec!wjb0414
        REAL,DIMENSION(2, 3000) :: GCD !����������������,�� �����꣬�̣߳� !wjb0414
    End type Canal
    !////////////////////////////////////////
    Type(Canal),ALLOCATABLE::Ca_cl(:)   !canal class�趨��������
    Type(Canal),ALLOCATABLE::Ca_cl_New(:)   !canal class�趨��������
    
    
    !/////���������////////////////
    Type::Areachange
        INTEGER::Serial_NUM !�ֲ����
        CHARACTER(len=12)::Name !�ڵ�����
        REAL::Incoord !������
        REAL::Outcoord !�������
        REAL::Locallosscoeff    !�ֲ�ˮ����ʧϵ��
        REAL::Additionlosscoeff  !������ʧϵ��
        REAL::SideSlopes      !����ϵ��
    End type Areachange
    !/////////////////////////////////////////
    Type(Areachange),ALLOCATABLE::Arch_cl(:)    !areachange class �趨���������
    
    !/////����������////////////////
    Type::Invertedsiphon
        INTEGER::Serial_NUM !�ֲ����
        CHARACTER(len=12)::Name !�ڵ�����
        REAL::Incoord !������
        REAL::Outcoord !�������
        INTEGER::ParallelNum    !���Ÿ���
        REAL::BottomWidth   !���
        REAL::Height    !�߶�
        REAL::Length    !����
        REAL::InInvertElev  !��ڵ׸߳�
        REAL::OutInvertElev  !���ڵ׸߳�
        REAL::Roughness !����������
        REAL::Insi_coeff !�ֲ�ˮ����ʧϵ��
    End type Invertedsiphon
    !//////////////////////////////////////////////
    Type(Invertedsiphon),ALLOCATABLE::Insi_cl(:)    !Invertedsiphon class �趨����������
    
    !////�Ŷ�����//////////
    Type::Bridge
        INTEGER::Serial_NUM !�ֲ����
        CHARACTER(len=12)::Name !�ڵ�����
        REAL::Incoord !������
        REAL::Outcoord !�������
        REAL::BottomWidth   !�׿�
        REAL::SideSlopes    !����ϵ��
        REAL::InInvertElev  !��ڵ׸߳�
        REAL::OutInvertElev  !���ڵ׸߳�
        REAL::Roughness !��������
        REAL::HinderWidth !����ˮ���
    End type Bridge
    !////////////////////////////////////////
    Type(Bridge),ALLOCATABLE::Brid_cl(:)   !Bridge class�趨��������
    
    !/////����բ����////////////////
    Type::Controlgate
        INTEGER::Serial_NUM !�ֲ����
        CHARACTER(len=12)::Name !�ڵ�����
        REAL::Incoord !������
        REAL::Outcoord !�������
        REAL::Initial_gate_Zup  !��ʼˮλ
        INTEGER::ParallelNum    !���Ÿ���
        REAL::SingleWidth   !���
        REAL::InvertElev  !�׸߳�
        REAL::Opendegree    !բ�ſ���
        INTEGER:: RUNNum    !��������
        REAL::G_A,G_B,G_C      !��բ����ϵ����բ�ſ���֮��Ĺ�ϵ����
        REAL::Discoeff      !բ�Ź���ϵ��
        REAL::Lenth
    End type Controlgate
    !//////////////////////////////////////////////
    Type(Controlgate),ALLOCATABLE::Cogate_cl(:)    !Controlgate class �趨����բ����
!    REAL,ALLOCATABLE::Opendegree(:,:)   !բ�ſ��ȹ���
!    REAL,ALLOCATABLE::C_Opendegree(:,:)   !�������ʹ��բ�ſ���

    !/////��ˮ������////////////////
    Type::SideOutFlow
        INTEGER::Serial_NUM !�ֲ����
        CHARACTER(len=12)::Name !�ڵ�����
        REAL::Incoord !������
        REAL::Outcoord !�������
        REAL::OutDischarge  !��ˮ����
    End type SideOutFlow
    !//////////////////////////////////////////////
    Type(SideOutFlow),ALLOCATABLE::SideOF_cl(:)    !SideOutFlow class �趨����բ����
    REAL,ALLOCATABLE::OutDischarge(:,:)     !��ˮ�ڳ�������
    REAL,ALLOCATABLE::C_OutDischarge(:,:)     !�������ʹ�÷�ˮ�ڳ���
    
    !/////��վ����////////////////
    Type::Pumpingstation
        INTEGER::Serial_NUM !�ֲ����
        CHARACTER(len=12)::Name !�ڵ�����
        REAL::Incoord !������
        REAL::Outcoord !�������
        REAL::Initial_pump_Zup  !��ʼˮλ
        REAL::InvertElev  !�׸߳�
        INTEGER::UNIT_NUM   !��վ�������
        INTEGER,ALLOCATABLE::Pump_Type(:)   !��վ��ˮ�����͵ı��
        REAL,ALLOCATABLE::Coefficient_A(:)   !������ϵ������
        REAL,ALLOCATABLE::Coefficient_B(:)    !һ����ϵ������
        REAL,ALLOCATABLE::Coefficient_C(:)    !����������
    End type Pumpingstation
    !//////////////////////////////////////////////
    Type(Pumpingstation),ALLOCATABLE::PumpST_cl(:)    !Pumpingstation class �趨��վ����
    Type(Pumpingstation),ALLOCATABLE::Ini_PumpST_cl(:)
    
    !/////ˮ������////////////////
    Type::Pump
        INTEGER::Serial_NUM !���ͱ��
        CHARACTER(len=12)::Name !ˮ������
        INTEGER::Coefficient_NUM   !ת�Ǹ���
        REAL,ALLOCATABLE::Bladeangle_Serial(:)   !ת������
        REAL,ALLOCATABLE::Coefficient_A(:)   !������ϵ������
        REAL,ALLOCATABLE::Coefficient_B(:)    !һ����ϵ������
        REAL,ALLOCATABLE::Coefficient_C(:)    !����������
    END TYPE Pump
    !//////////////////////////////////////////////
    Type(Pump),ALLOCATABLE::Pump_cl(:)    !Pumpingstation class �趨ˮ������
    
    !/////ת�ǹ���ʱ������//////////
    Type::Bladeangle_TimeSerial
        REAL,ALLOCATABLE::Each_PS(:,:)   !��վ��ÿ���õ�ת�ǹ���   ��ת�ǣ�ʱ�䣩
    END TYPE Bladeangle_TimeSerial
    Type(Bladeangle_TimeSerial),ALLOCATABLE::Bladeangle_cl(:)       !�������Ǳ�վ���
    Type(Bladeangle_TimeSerial),ALLOCATABLE::C_Bladeangle_cl(:)       !�������Ǳ�վ���
    
    !//////���ѭ������
    INTEGER::Cycle_NUM  !ÿһ�����������ѭ����
    !///////���Է������е�δ֪��
    
    !/////��������
    INTEGER::endflag=20  !5�������������������
    
    !/////��ʼ����
    INTEGER::startflag=1
    
    !/////����������
    INTEGER::bestoutputflag=1
    
    !////////�ɵ�����������ģ���ȶ��Բ����Ͳ��ʡ��ֲ�ˮ����ʧ�ȣ�
    REAL::Theta !�ռ�Ȩ��ϵ��(��������)
    REAL::Phi   !ʱ��Ȩ��ϵ��
    REAL::Psi   !�ռ�Ȩ��ϵ��(��������)
    REAL::PsiR  !�ռ�Ȩ��ϵ��(Ħ���¶�)
    REAL::Beta=1.0  !���涯��У��ϵ��
    REAL::GRAV=9.81 !�������ٶ�
    REAL::RoughnessCoeff    !���� ��������������ͬһ������ϵ����
    REAL,ALLOCATABLE::Insi_coeff(:,:)   !��������ʧϵ����Ŀǰ����2����ʧϵ����������Ҫ���е�����
    REAL::Discoeff_GATE !բ�Ź���ϵ��


    !///////���̸���
    INTEGER::Eq_NUM 
    !///////����ϵ������

    !////ͬ���۲���
    
    !///////////����ͬ�������˲�����ر���˵��
    INTEGER::NUMMember !�����г�Ա����

    INTEGER::Boundary_sign_UP,Boundary_sign_DOWN  !�߽��������Ʒ���0-������1-ˮλ��2-ˮλ~������ϵ
    INTEGER::Output_sign    !���������Ʒ���0-ֻ�������բ�������1-�����Ԫ��������
    REAL::m_coeff    !ĩ��բ�Ź�բ����ϵ������
    INTEGER::Control_sign  !ͬ��ģ�ͺ�ģ��ģ�͵Ŀ��Ʒ���ģ��Ϊ0��ͬ��Ϊ1
    INTEGER::Statetype   !״̬�����ͣ�0-���ʣ�1-����2-����+����
    INTEGER::Obsertype   !�۲������ͣ�0-ˮλ��1-������2-ˮλ+������
    REAL::PPQ, VVQ  !ˮλ~������ϵ����ϵ��
    


    !////////�Ż�����ʹ�ñ���
    !////////��վ
    REAL::Coefficient_1_5(34,4) !�洢��1-5����վ��ת��ϵ��,�Ѿ���������һά����ţ��ڶ�-��ά��ϵ�����ֱ��ӦA,B,C��
    INTEGER::Bladeangle_1_5(34,3) !ÿһ��ת��ϵ����Ӧ��ת�ǣ���һ-��ά�ֱ��ʾ1��2��3�Ż����ת��ֵ
    REAL::Coefficient_6(55,4)   !�洢��6����վ��ת��ϵ�����Ѿ���������һά����ţ��ڶ�-��ά��ϵ�����ֱ��ӦA,B,C��
    INTEGER::Bladeangle_6(55,3)     !ÿһ��ת��ϵ����Ӧ��ת�ǣ���һ-��ά�ֱ��ʾ1��2��3�Ż����ת��ֵ
    REAL,ALLOCATABLE::Forebay_MAX(:),Forebay_MIN(:),Outlet_MAX(:),Outlet_MIN(:),DESIGN_ZIN(:) !�洢������վվǰ��վ����ߡ����ˮλԼ��|   
    INTEGER::Bladeangle_varibale(30)

    !����բ
!    REAL,ALLOCATABLE::ST_STATE(:,:)     !�洢������բÿСʱբǰˮλ��բ��ˮλ
!    REAL,ALLOCATABLE::DS_STATE(:)
!    INTEGER,ALLOCATABLE::SG_POSTION(:)	!����բλ��
!     INTEGER,ALLOCATABLE::DM_POSTION(:)  !����բ+��β����λ��
!     INTEGER,ALLOCATABLE::ELE_POSTION(:) !�������Ŷ������Ԫ��λ��
!    REAL,ALLOCATABLE::GAV_Opendegree(:)
!    REAL,ALLOCATABLE::Deltadegree(:,:)
    REAL,ALLOCATABLE::ZQSW_TARGET(:)    !����բբǰĿ��ˮλ���м�������բ���м���
    REAL::D_WL_TARGET   !�����α߽紦��Ŀ��ˮλ
    REAL,ALLOCATABLE::ST_CON_UP(:)  !ˮλԼ������
    REAL,ALLOCATABLE::ST_CON_DOWN(:)    !ˮλԼ������
    REAL::DS_CON_UP,DS_CON_DOWN !����ˮλ�߽磬���ޡ�����
    REAL,ALLOCATABLE::h_new_INIT(:),QQ_INIT(:,:)    !��ʼˮλ����ʼ����
    REAL,ALLOCATABLE::Opendegree_INIT(:,:)  !��ʼ����
    
    
    REAL,ALLOCATABLE::BestOpendegree(:,:)
    REAL,ALLOCATABLE::BestST_STATE(:,:)
    REAL,ALLOCATABLE::BestQControlGate(:,:)
    REAL,ALLOCATABLE::BestVolume(:,:)
End


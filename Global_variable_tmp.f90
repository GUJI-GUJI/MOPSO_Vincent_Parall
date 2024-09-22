MODULE GLOBAL_VARIABLE_TMP  !�������
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    !����ȫ�ֱ���
    !�����ļ�
    INTEGER::Start_Time,End_Time   !��ʼʱ�䣬����ʱ�䣨�����ʾ��
    INTEGER::Boundary_Step,Boundary_NUM !�߽������������߽����и������߽������еĹ۲����ݣ�
    INTEGER::DeltaT,Calculate_NUM,Out_DeltaT !����ʱ�䲽����������������ʱ�䲽��
    REAL::DeltaX    !�ռ䲽��
    REAL,ALLOCATABLE::delX(:) !�洢ÿһ�������ķֶ���
    INTEGER::c_element    !��������ռ䲽����������Ԫ���ĸ���
    INTEGER,ALLOCATABLE::n(:)   !ÿһ��Ԫ���ķֶ����
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

    !/////�����Ϣ
    REAL,ALLOCATABLE::h(:),Q(:),Z(:),vc(:),L(:),Vol(:),Pool_Volume_ini(:)
    !////////////////////////////////////////
    Type(Canal),ALLOCATABLE::Ca_cl(:)   !canal class�趨��������
    Type(Canal),ALLOCATABLE::Ca_cl_New(:)   !canal class�趨��������
    
 
    !/////////////////////////////////////////
    Type(Areachange),ALLOCATABLE::Arch_cl(:)    !areachange class �趨���������
    

    !//////////////////////////////////////////////
    Type(Invertedsiphon),ALLOCATABLE::Insi_cl(:)    !Invertedsiphon class �趨����������
    

    !////////////////////////////////////////
    Type(Bridge),ALLOCATABLE::Brid_cl(:)   !Bridge class�趨��������
    

    !//////////////////////////////////////////////
    Type(Controlgate),ALLOCATABLE::Cogate_cl(:)    !Controlgate class �趨����բ����
    REAL,ALLOCATABLE::Opendegree(:,:)   !բ�ſ��ȹ���
    REAL,ALLOCATABLE::C_Opendegree(:,:)   !�������ʹ��բ�ſ���
    REAL,ALLOCATABLE::C_Bladeangle(:,:,:)   !�������ʹ��բ�ſ���

 
    !//////////////////////////////////////////////
    Type(SideOutFlow),ALLOCATABLE::SideOF_cl(:)    !SideOutFlow class �趨����բ����
    REAL,ALLOCATABLE::OutDischarge(:,:)     !��ˮ�ڳ�������
    REAL,ALLOCATABLE::C_OutDischarge(:,:)     !�������ʹ�÷�ˮ�ڳ���
    

    !//////////////////////////////////////////////
    Type(Pumpingstation),ALLOCATABLE::PumpST_cl(:)    !Pumpingstation class �趨��վ����
    Type(Pumpingstation),ALLOCATABLE::Ini_PumpST_cl(:)
    

    !//////////////////////////////////////////////
    Type(Pump),ALLOCATABLE::Pump_cl(:)    !Pumpingstation class �趨ˮ������
    

    Type(Bladeangle_TimeSerial),ALLOCATABLE::Bladeangle_cl(:)       !�������Ǳ�վ���
    Type(Bladeangle_TimeSerial),ALLOCATABLE::C_Bladeangle_cl(:)       !�������Ǳ�վ���
    
    !//////���ѭ������
    INTEGER::Cycle_NUM  !ÿһ�����������ѭ����
    !///////���Է������е�δ֪��
    
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

    !/////////////�Ǻ㶨������
    !δ֪���洢
    REAL,ALLOCATABLE::ZZ(:,:)
    REAL,ALLOCATABLE::QQ(:,:)
    REAL,ALLOCATABLE::Pool_Volume(:,:)
    REAL,ALLOCATABLE::Total_Volume(:)
    !��һ����һֱ���洢
    REAL,ALLOCATABLE::ZZOLD(:,:)
    REAL,ALLOCATABLE::QQOLD(:,:)
    REAL,ALLOCATABLE::ELE_Volume(:,:)
    !׷��ϵ��
    REAL,ALLOCATABLE::SS(:)
    REAL,ALLOCATABLE::TT(:)
    REAL,ALLOCATABLE::PP(:)
    REAL,ALLOCATABLE::VV(:)
    !��̬������ˮ�����飨�£�
    REAL,ALLOCATABLE::h_new(:)
    !///////���̸���
    INTEGER::Eq_NUM 
    !///////����ϵ������

    !////ͬ���۲���
    !REAL,ALLOCATABLE::Si_z(:,:)   !ˮλģ��ֵ����
    !REAL,ALLOCATABLE::Si_q(:,:)   !����ģ��ֵ����
    !REAL,ALLOCATABLE::Si_v(:,:)   !����ģ��ֵ����
    

    !///////////����ͬ�������˲�����ر���˵��
    INTEGER::NUMMember !�����г�Ա����

    INTEGER::Boundary_sign_UP,Boundary_sign_DOWN  !�߽��������Ʒ���0-������1-ˮλ��2-ˮλ~������ϵ
    INTEGER::Output_sign    !���������Ʒ���0-ֻ�������բ�������1-�����Ԫ��������
    REAL::m_coeff    !ĩ��բ�Ź�բ����ϵ������
    INTEGER::Control_sign  !ͬ��ģ�ͺ�ģ��ģ�͵Ŀ��Ʒ���ģ��Ϊ0��ͬ��Ϊ1
    INTEGER::Statetype   !״̬�����ͣ�0-���ʣ�1-����2-����+����
    INTEGER::Obsertype   !�۲������ͣ�0-ˮλ��1-������2-ˮλ+������
    REAL::PPQ, VVQ  !ˮλ~������ϵ����ϵ��
    
    INTEGER,ALLOCATABLE::SG_POSTION(:)  !����բλ��
    INTEGER,ALLOCATABLE::DM_POSTION(:),DDM_POSTION(:)  !����բ+��β����λ��
    INTEGER,ALLOCATABLE::ELE_POSTION(:) !�������Ŷ������Ԫ��λ��
End


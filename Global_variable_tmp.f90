MODULE GLOBAL_VARIABLE_TMP  !定义变量
    USE GLOBAL_VARIABLE
    IMPLICIT NONE
    !定义全局变量
    !输入文件
    INTEGER::Start_Time,End_Time   !初始时间，结束时间（以秒表示）
    INTEGER::Boundary_Step,Boundary_NUM !边界条件步长，边界序列个数（边界条件中的观测数据）
    INTEGER::DeltaT,Calculate_NUM,Out_DeltaT !计算时间步长，计算次数，输出时间步长
    REAL::DeltaX    !空间步长
    REAL,ALLOCATABLE::delX(:) !存储每一个渠道的分段数
    INTEGER::c_element    !计算调整空间步长后参与计算元件的个数
    INTEGER,ALLOCATABLE::n(:)   !每一个元件的分段情况
    INTEGER,ALLOCATABLE::n_canal(:) !每个渠道的分段情况
    REAL,ALLOCATABLE::R_Upstream(:),R_Downstream(:) !上游边界条件和下游边界条件
    REAL,ALLOCATABLE::Upstream(:),Downstream(:) !上游边界条件和下游边界条件
    INTEGER::Element_NUM    !计算元件的个数
    INTEGER::Canal_NUM,Areachange_NUM,Invertedsiphon_NUM,Controlgate_NUM,SideOutFlow_NUM,Bridge_NUM,PumpST_NUM    !渠道个数，渐变段个数，倒虹吸个数，闸门个数,分水口个数,桥墩个数，泵站个数
    INTEGER::Pump_type_NUM !水泵类型个数
    INTEGER::Canal_New_NUM  !实际参与计算渠段的个数
    !用系统编号标识计算元件在整个河网中的顺序。局部编号为在同类元件中的序号。类型代码表示计算元件的类型
    INTEGER,ALLOCATABLE::Ele_Rel(:,:)  !element_relation存储系统编号、类型代码、局部编号
    INTEGER,ALLOCATABLE::Ele_Rel_New(:,:)   !element_relation_new存储考虑时间步长和空间步长之后的系统编号、类型代码、局部编号
    REAL:: Initial_Discharge    !上游边界也选择水位时的初始流量
    REAL:: Initial_gate_Zup     !初始闸前水位
    REAL:: Initial_pump_Zup     !初始站前水位
    REAL:: D_elevation  !下游边界底高程
    REAL::Waterdepth(10000,2),Waterdepth_flag
    REAL,ALLOCATABLE::Dynamic_Result_Zup(:,:),Dynamic_Result_Zdown(:,:),Dynamic_Result_Q(:,:),Dynamic_Result_Zbd(:,:),Dynamic_Result_Qbd(:,:),Dynamic_Result_ZZ(:,:),Dynamic_Result_QQ(:,:)

    !/////结果信息
    REAL,ALLOCATABLE::h(:),Q(:),Z(:),vc(:),L(:),Vol(:),Pool_Volume_ini(:)
    !////////////////////////////////////////
    Type(Canal),ALLOCATABLE::Ca_cl(:)   !canal class设定渠道类型
    Type(Canal),ALLOCATABLE::Ca_cl_New(:)   !canal class设定渠道类型
    
 
    !/////////////////////////////////////////
    Type(Areachange),ALLOCATABLE::Arch_cl(:)    !areachange class 设定渐变段类型
    

    !//////////////////////////////////////////////
    Type(Invertedsiphon),ALLOCATABLE::Insi_cl(:)    !Invertedsiphon class 设定倒虹吸类型
    

    !////////////////////////////////////////
    Type(Bridge),ALLOCATABLE::Brid_cl(:)   !Bridge class设定渠道类型
    

    !//////////////////////////////////////////////
    Type(Controlgate),ALLOCATABLE::Cogate_cl(:)    !Controlgate class 设定节制闸类型
    REAL,ALLOCATABLE::Opendegree(:,:)   !闸门开度过程
    REAL,ALLOCATABLE::C_Opendegree(:,:)   !计算过程使用闸门开度
    REAL,ALLOCATABLE::C_Bladeangle(:,:,:)   !计算过程使用闸门开度

 
    !//////////////////////////////////////////////
    Type(SideOutFlow),ALLOCATABLE::SideOF_cl(:)    !SideOutFlow class 设定节制闸类型
    REAL,ALLOCATABLE::OutDischarge(:,:)     !分水口出流过程
    REAL,ALLOCATABLE::C_OutDischarge(:,:)     !计算过程使用分水口出流
    

    !//////////////////////////////////////////////
    Type(Pumpingstation),ALLOCATABLE::PumpST_cl(:)    !Pumpingstation class 设定泵站类型
    Type(Pumpingstation),ALLOCATABLE::Ini_PumpST_cl(:)
    

    !//////////////////////////////////////////////
    Type(Pump),ALLOCATABLE::Pump_cl(:)    !Pumpingstation class 设定水泵类型
    

    Type(Bladeangle_TimeSerial),ALLOCATABLE::Bladeangle_cl(:)       !括号里是泵站编号
    Type(Bladeangle_TimeSerial),ALLOCATABLE::C_Bladeangle_cl(:)       !括号里是泵站编号
    
    !//////最大循环次数
    INTEGER::Cycle_NUM  !每一迭代步的最大循环数
    !///////线性方程组中的未知量
    
    !////////可调参数（包括模型稳定性参数和糙率、局部水力损失等）
    REAL::Theta !空间权重系数(连续方程)
    REAL::Phi   !时间权重系数
    REAL::Psi   !空间权重系数(动量方程)
    REAL::PsiR  !空间权重系数(摩擦坡度)
    REAL::Beta=1.0  !断面动量校正系数
    REAL::GRAV=9.81 !重力加速度
    REAL::RoughnessCoeff    !糙率 （沿线渠道采用同一个糙率系数）
    REAL,ALLOCATABLE::Insi_coeff(:,:)   !倒虹吸损失系数（目前设置2个损失系数，后期需要进行调整）
    REAL::Discoeff_GATE !闸门过流系数

    !/////////////非恒定流计算
    !未知量存储
    REAL,ALLOCATABLE::ZZ(:,:)
    REAL,ALLOCATABLE::QQ(:,:)
    REAL,ALLOCATABLE::Pool_Volume(:,:)
    REAL,ALLOCATABLE::Total_Volume(:)
    !上一步的一直量存储
    REAL,ALLOCATABLE::ZZOLD(:,:)
    REAL,ALLOCATABLE::QQOLD(:,:)
    REAL,ALLOCATABLE::ELE_Volume(:,:)
    !追赶系数
    REAL,ALLOCATABLE::SS(:)
    REAL,ALLOCATABLE::TT(:)
    REAL,ALLOCATABLE::PP(:)
    REAL,ALLOCATABLE::VV(:)
    !稳态各断面水深数组（新）
    REAL,ALLOCATABLE::h_new(:)
    !///////方程个数
    INTEGER::Eq_NUM 
    !///////定义系数矩阵

    !////同化观测量
    !REAL,ALLOCATABLE::Si_z(:,:)   !水位模拟值序列
    !REAL,ALLOCATABLE::Si_q(:,:)   !流量模拟值序列
    !REAL,ALLOCATABLE::Si_v(:,:)   !流速模拟值序列
    

    !///////////数据同化集合滤波的相关变量说明
    INTEGER::NUMMember !集合中成员个数

    INTEGER::Boundary_sign_UP,Boundary_sign_DOWN  !边界条件控制符，0-流量，1-水位，2-水位~流量关系
    INTEGER::Output_sign    !输出结果控制符，0-只输出节制闸处结果；1-输出各元件断面结果
    REAL::m_coeff    !末端闸门过闸流量系数计算
    INTEGER::Control_sign  !同化模型和模拟模型的控制符，模拟为0，同化为1
    INTEGER::Statetype   !状态量类型（0-糙率；1-局损；2-糙率+局损）
    INTEGER::Obsertype   !观测量类型（0-水位；1-流量；2-水位+流量）
    REAL::PPQ, VVQ  !水位~流量关系断面系数
    
    INTEGER,ALLOCATABLE::SG_POSTION(:)  !节制闸位置
    INTEGER,ALLOCATABLE::DM_POSTION(:),DDM_POSTION(:)  !节制闸+首尾断面位置
    INTEGER,ALLOCATABLE::ELE_POSTION(:) !渠道、桥墩以外的元件位置
End


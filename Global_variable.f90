    MODULE GLOBAL_VARIABLE  !定义变量
    IMPLICIT NONE
    !定义全局变量

    !V2.0V2.0V2.0
    !V2.0V2.0V2.0
    !V2.0V2.0V2.0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!基本常数!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    REAL,PARAMETER :: PI=3.14159265

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!工程参数类型!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !////渠道类型//////////
    Type::Canal
        INTEGER::Serial_NUM                 !局部编号
        INTEGER::irr_Serial_NUM=0           !表示该渠道是不规则渠道中的第几条
        CHARACTER(len=12)::Name             !节点名称
        REAL::Incoord                       !入口里程
        REAL::Outcoord                      !出口里程
        REAL::BottomWidth                   !底宽
        REAL::SideSlopes                    !边坡系数
        REAL::InInvertElev                  !入口底高程
        REAL::OutInvertElev                 !出口底高程
        REAL::InTopElev                     !入口顶高程
        REAL::OutTopElev                    !出口顶高程
        REAL::Roughness                     !渠道糙率
        INTEGER :: CroSecTYPE               !断面类型，0 梯形（矩形），1圆形，2不规则断面
        INTEGER  :: CDS_in                  !#TU入口断面测点数
        INTEGER  :: CDS_out                 !#TU入口断面测点数
        REAL,DIMENSION(2, 3000) :: GCD_in   !#TU入口高程对
        REAL,DIMENSION(2, 3000) :: GCD_out  !#TU出口高程对
    End type Canal

    !/////渐变段类型////////////////
    Type::Areachange
        INTEGER::Serial_NUM         !局部编号
        CHARACTER(len=12)::Name     !节点名称
        REAL::Incoord               !入口里程
        REAL::Outcoord              !出口里程
        REAL::Locallosscoeff        !局部水力损失系数
        REAL::Additionlosscoeff     !其他损失系数
        REAL::SideSlopes            !边坡系数
    End type Areachange

    !/////倒虹吸类型////////////////
    Type::Invertedsiphon
        INTEGER::Serial_NUM     !局部编号
        CHARACTER(len=12)::Name !节点名称
        REAL::Incoord           !入口里程
        REAL::Outcoord          !出口里程
        INTEGER::ParallelNum    !并排个数
        REAL::BottomWidth       !宽度
        REAL::Height            !高度
        REAL::Length            !长度
        REAL::InInvertElev      !入口底高程
        REAL::OutInvertElev     !出口底高程
        REAL::Roughness         !倒虹吸糙率
        REAL::Insi_coeff        !局部水力损失系数
    End type Invertedsiphon

    !////桥墩类型//////////
    Type::Bridge
        INTEGER::Serial_NUM     !局部编号
        CHARACTER(len=12)::Name !节点名称
        REAL::Incoord           !入口里程
        REAL::Outcoord          !出口里程
        REAL::BottomWidth       !底宽
        REAL::SideSlopes        !边坡系数
        REAL::InInvertElev      !入口底高程
        REAL::OutInvertElev     !出口底高程
        REAL::Roughness         !渠道糙率
        REAL::HinderWidth       !总阻水宽度
    End type Bridge

    !/////节制闸类型////////////////
    Type::Controlgate
        INTEGER::Serial_NUM     !局部编号
        CHARACTER(len=12)::Name !节点名称
        REAL::Incoord           !入口里程
        REAL::Outcoord          !出口里程
        REAL::Initial_gate_Zup  !初始水位
        INTEGER::ParallelNum    !并排个数
        REAL::SingleWidth       !宽度
        REAL::InvertElev        !底高程
        REAL::Opendegree        !闸门开度
        INTEGER:: RUNNum        !开启孔数
        REAL::G_A,G_B,G_C       !过闸流量系数与闸门开度之间的关系曲线
        REAL::Discoeff          !闸门过流系数
        REAL::Lenth
    End type Controlgate

    !/////分水口类型////////////////
    Type::SideOutFlow
        INTEGER::Serial_NUM     !局部编号
        CHARACTER(len=12)::Name !节点名称
        REAL::Incoord           !入口里程
        REAL::Outcoord          !出口里程
        REAL::OutDischarge      !分水流量
    End type SideOutFlow

    !/////泵站类型////////////////
    Type::Pumpingstation
        INTEGER::Serial_NUM                 !局部编号
        CHARACTER(len=12)::Name             !节点名称
        REAL::Incoord                       !入口里程
        REAL::Outcoord                      !出口里程
        REAL::Initial_pump_Zup              !初始水位
        REAL::InvertElev                    !底高程
        INTEGER::UNIT_NUM                   !泵站机组个数
        INTEGER,ALLOCATABLE::Pump_Type(:)   !泵站中水泵类型的编号
        REAL,ALLOCATABLE::Coefficient_A(:)  !二次项系数序列
        REAL,ALLOCATABLE::Coefficient_B(:)  !一次项系数序列
        REAL,ALLOCATABLE::Coefficient_C(:)  !常数项序列
    End type Pumpingstation

    !/////水泵类型////////////////
    Type::Pump
        INTEGER::Serial_NUM                     !类型编号
        CHARACTER(len=12)::Name                 !水泵名称
        INTEGER::Coefficient_NUM                !转角个数
        REAL,ALLOCATABLE::Bladeangle_Serial(:)  !转角序列
        REAL,ALLOCATABLE::Coefficient_A(:)      !二次项系数序列
        REAL,ALLOCATABLE::Coefficient_B(:)      !一次项系数序列
        REAL,ALLOCATABLE::Coefficient_C(:)      !常数项序列
    END TYPE Pump

    !/////转角过程时间序列//////////
    Type::Bladeangle_TimeSerial
        REAL,ALLOCATABLE::Each_PS(:,:)   !泵站中每个泵的转角过程   （转角，时间）
    END TYPE Bladeangle_TimeSerial

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!工况参数!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !边界条件
    TYPE::Boundary
        CHARACTER(LEN=50)::name             !名称
        CHARACTER(LEN=50)::Filename         !文件地址
        INTEGER::Start_Time                 !初始时间
        INTEGER::End_Time                   !结束时间（以秒表示）
        INTEGER::Step                       !边界条件步长
        INTEGER::Upsign                     !边界条件控制符，0-流量，1-水位，2-水位~流量关系
        INTEGER::Downsign                   !边界条件控制符，0-流量，1-水位，2-水位~流量关系
        REAL,ALLOCATABLE::UpQ(:)            !上游流量过程
        REAL,ALLOCATABLE::UpZ(:)            !上游水位过程
        REAL,ALLOCATABLE::DownQ(:)          !下游流量过程
        REAL,ALLOCATABLE::DownZ(:)          !下游水位过程
        REAL,ALLOCATABLE::Opendegree(:,:)   !闸门开度过程
        REAL,ALLOCATABLE::OutDischarge(:,:) !分水口出流过程
        Type(Bladeangle_TimeSerial),&
            ALLOCATABLE::Bladeangle_cl(:)   !括号里是泵站编号

        LOGICAL::Auto_Update                !是否自动更新下列参数（若自动更新则在check中进行更新）
        INTEGER::NUM                        !边界序列个数（边界条件中的观测数据）
        INTEGER::Gate_NUM                   !闸门数量
        INTEGER::PumpST_NUM                !泵站数量
        INTEGER::SideOutFlow_NUM            !分水口数量
    END TYPE Boundary

    !工程信息
    TYPE::Com_element
        CHARACTER(LEN=50)::name                         !名称
        CHARACTER(LEN=50)::Filename                     !文件地址
        INTEGER::Element_NUM                            !计算元件的个数
        INTEGER,ALLOCATABLE::Ele_Rel(:,:)               !element_relation存储系统编号、类型代码、局部编号
        INTEGER::Canal_NUM                              !渠道个数
        Type(Canal),ALLOCATABLE::Ca_cl(:)               !canal class设定渠道类型
        INTEGER::Areachange_NUM                         !渐变段个数
        Type(Areachange),ALLOCATABLE::Arch_cl(:)        !areachange class 设定渐变段类型
        INTEGER::Invertedsiphon_NUM                     !倒虹吸个数
        Type(Invertedsiphon),ALLOCATABLE::Insi_cl(:)    !Invertedsiphon class 设定倒虹吸类型
        INTEGER::Controlgate_NUM                        !闸门个数
        Type(Controlgate),ALLOCATABLE::Cogate_cl(:)     !Controlgate class 设定节制闸类型
        INTEGER::SideOutFlow_NUM                        !分水口个数
        Type(SideOutFlow),ALLOCATABLE::SideOF_cl(:)     !SideOutFlow class 设定分水口类型
        INTEGER::Bridge_NUM                             !桥墩个数
        Type(Bridge),ALLOCATABLE::Brid_cl(:)            !Bridge class设定渠道类型
        INTEGER::PumpST_NUM                             !泵站个数
        Type(Pumpingstation),ALLOCATABLE::PumpST_cl(:)  !Pumpingstation class 设定泵站类型
        INTEGER::Pump_type_NUM                          !水泵类型个数
        Type(Pump),ALLOCATABLE::Pump_cl(:)              !Pumpingstation class 设定水泵类型

        LOGICAL::Auto_Update                            !是否自动更新下列参数（若自动更新则在check中进行更新）
        INTEGER,ALLOCATABLE::n(:)                       !每一个元件的分段情况
        INTEGER,ALLOCATABLE::n_canal(:)                 !每个渠道的分段情况
        INTEGER,ALLOCATABLE::SG_POSTION(:)              !节制闸位置
        INTEGER,ALLOCATABLE::DM_POSTION(:)              !节制闸+首尾断面位置
        INTEGER,ALLOCATABLE::ELE_POSTION(:)             !渠道、桥墩以外的元件位置
    END TYPE Com_element

    !对应时刻的各参数值
    TYPE::Condition
        CHARACTER(LEN=50)::Comname          !对应工程信息名
        CHARACTER(LEN=50)::Boundaryname     !对应边界条件信息名
        REAL::time                          !时刻
        REAL::h(10000)                      !初始水深，原waterdepth
        REAL::Q(10000)                      !初始流量，原waterdepth
        REAL,ALLOCATABLE::ZZ(:,:)           !水位
        REAL,ALLOCATABLE::QQ(:,:)           !流量
        REAL,ALLOCATABLE::ELE_Volume(:,:)   !蓄量
        REAL,ALLOCATABLE::SS(:)             !追赶系数
        REAL,ALLOCATABLE::TT(:)             !追赶系数
        REAL,ALLOCATABLE::PP(:)             !追赶系数
        REAL,ALLOCATABLE::VV(:)             !追赶系数

        LOGICAL::Auto_Update                !是否自动计算下列参数
        REAL,ALLOCATABLE::Pool_Volume(:,:)  !渠池蓄量
        REAL,ALLOCATABLE::Total_Volume(:)   !总蓄量
    END TYPE Condition

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!模型参数类型!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !基本计算参数，包括模型稳定性参数和糙率、局部水力损失等
    TYPE::Coeff
        REAL::Theta                         !空间权重系数(连续方程)
        REAL::Phi                           !时间权重系数
        REAL::Psi                           !空间权重系数(动量方程)
        REAL::PsiR                          !空间权重系数(摩擦坡度)
        REAL::Beta                          !断面动量校正系数
        REAL::GRAV                          !重力加速度
        REAL::RoughnessCoeff                !糙率 （沿线渠道采用同一个糙率系数）
        REAL,ALLOCATABLE::Insi_coeff(:,:)   !倒虹吸损失系数（目前设置2个损失系数，后期需要进行调整）
        REAL::Discoeff_GATE                 !闸门过流系数
    END TYPE Coeff

    !同化参数
    TYPE::Assimilate
        INTEGER::NUMMember      !集合中成员个数
        INTEGER::Output_sign    !输出结果控制符，0-只输出节制闸处结果；1-输出各元件断面结果
        REAL::m_coeff           !末端闸门过闸流量系数计算
        INTEGER::Control_sign   !同化模型和模拟模型的控制符，模拟为0，同化为1
        INTEGER::Statetype      !状态量类型（0-糙率；1-局损；2-糙率+局损）
        INTEGER::Obsertype      !观测量类型（0-水位；1-流量；2-水位+流量）
        REAL::PPQ               !水位~流量关系断面系数
        REAL::VVQ               !水位~流量关系断面系数
    END TYPE Assimilate

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!模型参数!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !模型总配置
    TYPE::Config
        CHARACTER(LEN=50)::Filename         !文件地址
        TYPE(Coeff)::Coeff_set
        TYPE(Assimilate)::Assimilate_set
    END TYPE Config

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!计算结果!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    !初始恒定流计算结果
    TYPE::Initial_Result
        CHARACTER(LEN=50)::Filename         !文件地址
        CHARACTER(LEN=50)::Comname          !对应工程信息名
        CHARACTER(LEN=50)::Boundaryname     !对应边界条件信息名
        REAL,ALLOCATABLE::h(:)              !水深
        REAL,ALLOCATABLE::Q(:)              !流量
        REAL,ALLOCATABLE::Z(:)              !水位
        REAL,ALLOCATABLE::vc(:)             !流速
        REAL,ALLOCATABLE::L(:)              !里程
        REAL,ALLOCATABLE::Vol(:)            !蓄量
        REAL,ALLOCATABLE::Pool_Volume(:)    !渠段蓄量

        LOGICAL::Auto_Update                !是否自动计算下列参数
        INTEGER::NUM                        !元件数量
    END TYPE Initial_Result

    !非恒定流结果
    TYPE::Dynamic_Result
        CHARACTER(LEN=50)::Filename         !文件地址
        CHARACTER(LEN=50)::Comname          !对应工程信息名
        CHARACTER(LEN=50)::Boundaryname     !对应边界条件信息名
        INTEGER::DeltaT                     !时间步长
        REAL,ALLOCATABLE::Zup(:,:)          !闸前水位
        REAL,ALLOCATABLE::Zdown(:,:)        !闸后水位
        REAL,ALLOCATABLE::Q(:,:)            !过闸流量
        REAL,ALLOCATABLE::Zbd(:,:)          !末端水位
        REAL,ALLOCATABLE::Qbd(:,:)          !末端流量
        REAL,ALLOCATABLE::ZZ(:,:)           !各断面水位
        REAL,ALLOCATABLE::QQ(:,:)           !各断面流量
        REAL,ALLOCATABLE::Pool_Volume(:,:)    !渠段蓄量
        REAL,ALLOCATABLE::Total_Volume(:)   !总蓄量

        INTEGER::NUM
    END TYPE Dynamic_Result

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    INTEGER::Out_DeltaT !计算时间步长，计算次数，输出时间步长
    INTEGER::Output_sign    !输出结果控制符，0-只输出节制闸处结果；1-输出各元件断面结果

    INTEGER::PumpST_NUM    !渠道个数，渐变段个数，倒虹吸个数，闸门个数,分水口个数,桥墩个数，泵站个数
    INTEGER::Pump_type_NUM !水泵类型个数
    !//////////////////////////////////////////////
    Type(Pumpingstation),ALLOCATABLE::PumpST_cl(:)    !Pumpingstation class 设定泵站类型
    Type(Pumpingstation),ALLOCATABLE::Ini_PumpST_cl(:)
    !//////////////////////////////////////////////
    Type(Pump),ALLOCATABLE::Pump_cl(:)    !Pumpingstation class 设定水泵类型

    REAL,ALLOCATABLE::C_Bladeangle(:,:,:)   !计算过程使用闸门开度
    Type(Bladeangle_TimeSerial),ALLOCATABLE::Bladeangle_cl(:)       !括号里是泵站编号
    Type(Bladeangle_TimeSerial),ALLOCATABLE::C_Bladeangle_cl(:)       !括号里是泵站编号






    !New——MPC
    REAL::MIN_Target1=10000,MIN_Target2=10000
    REAL,ALLOCATABLE::ZQSW_TARGET(:)    !节制闸闸前目标水位，有几个节制闸就有几个
    REAL::D_WL_TARGET   !最下游边界处的目标水位
    REAL,ALLOCATABLE::ST_CON_UP(:)  !水位约束上限
    REAL,ALLOCATABLE::ST_CON_DOWN(:)    !水位约束下限
    REAL::DS_CON_UP,DS_CON_DOWN !下游水位边界，上限、下限
    REAL,ALLOCATABLE::h_new_INIT(:),QQ_INIT(:,:)    !初始水位、初始流量
    !REAL,ALLOCATABLE::Opendegree_INIT(:,:)  !初始开度


    REAL,ALLOCATABLE::BestOpendegree(:,:)
    REAL,ALLOCATABLE::BestST_STATE(:,:)
    REAL,ALLOCATABLE::BestQControlGate(:,:)
    REAL,ALLOCATABLE::BestVolume(:,:)




    !//////最大循环次数
    INTEGER::Cycle_NUM  !每一迭代步的最大循环数
    !///////线性方程组中的未知量

    !/////收敛条件
    INTEGER::endflag=20  !5代连续满足条件则结束

    !/////开始条件
    INTEGER::startflag=1

    !/////最佳输出条件
    INTEGER::bestoutputflag=1





    !temp，这样用好像不好
    integer::Boudarys_Number
    integer::Controlgate_NUMber
    real,allocatable::Opendegree_INIT(:)!读取开度的第一个
    TYPE(Boundary)::Boundary_in,Boundary_cal,Boundary_2in,Boundary_2cal
    TYPE(Config)::Config_main
    TYPE(Com_element)::Com_element_main,Com_element_cal
    TYPE(Initial_Result)::Initial_Result_main
    TYPE(Condition)::Condition_0
    
    Integer::Function_index=0!evaluate函数的运行次数
    Real,allocatable::Ini_ArriveTime(:)!预测的污染到达时间
    Real::Ini_CloseTime!预测的闸门关闭最短所需时间
    Integer::total_BoundaryNum=97
    
    Integer::firstgate_Index
    Integer::firstGateIndex
    Integer::outflow_index
    Integer::arrivetime_weight
    Real::pollut_location!污染源距离上游距离
    Real::pollut_weight
    Real::pollut_ConveyDistance
    Real::outflow_maxQ
    Real::outflow_changelimit
    Real::inGate_changelimit
    INTEGER::generations
    INTEGER::populations
    !Integer::maxArchive,gens,pops,NP,NF
    End


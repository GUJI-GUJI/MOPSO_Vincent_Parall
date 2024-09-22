MODULE GLOBAL_VARIABLE  !定义变量
    IMPLICIT NONE
    !定义全局变量
    !输入文件
    INTEGER::Start_Time,End_Time   !初始时间，结束时间（以秒表示）
    INTEGER::Boundary_Step,Boundary_NUM !边界条件步长，边界序列个数（边界条件中的观测数据）
    INTEGER::Controltime    !调控次数
    INTEGER::DeltaT,Calculate_NUM,Out_DeltaT !计算时间步长，计算次数，输出时间步长
    REAL::DeltaX    !空间步长
    REAL,ALLOCATABLE::delX(:) !存储每一个渠道的分段数
    INTEGER::c_element    !计算调整空间步长后参与计算元件的个数
    INTEGER,ALLOCATABLE::nnn(:)   !每一个元件的分段情况
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
    REAL::MIN_Target1=10000,MIN_Target2=10000
    
    !/////结果信息
    REAL,ALLOCATABLE::h(:),Q(:),Z(:),vc(:),L(:),Vol(:)
    
    !////渠道类型//////////
    Type::Canal
        INTEGER::Serial_NUM !局部编号
        CHARACTER(len=12)::Name !节点名称
        REAL::Incoord !入口里程
        REAL::Outcoord !出口里程
        REAL::BottomWidth   !底宽
        REAL::SideSlopes    !边坡系数
        REAL::InInvertElev  !入口底高程
        REAL::OutInvertElev  !出口底高程
        REAL::InTopElev !入口顶高程
        REAL::OutTopElev    !出口顶高程
        REAL::Roughness !渠道糙率
        INTEGER :: CroSecTYPE  !断面类型，0 梯形（矩形），1圆形，2不规则断面 !wjb0414
        INTEGER  :: CDS  !断面测点数, 不规则断面编号 ID_IR_CroSec!wjb0414
        REAL,DIMENSION(2, 3000) :: GCD !不规则断面地形数组,（ 横坐标，高程） !wjb0414
    End type Canal
    !////////////////////////////////////////
    Type(Canal),ALLOCATABLE::Ca_cl(:)   !canal class设定渠道类型
    Type(Canal),ALLOCATABLE::Ca_cl_New(:)   !canal class设定渠道类型
    
    
    !/////渐变段类型////////////////
    Type::Areachange
        INTEGER::Serial_NUM !局部编号
        CHARACTER(len=12)::Name !节点名称
        REAL::Incoord !入口里程
        REAL::Outcoord !出口里程
        REAL::Locallosscoeff    !局部水力损失系数
        REAL::Additionlosscoeff  !其他损失系数
        REAL::SideSlopes      !边坡系数
    End type Areachange
    !/////////////////////////////////////////
    Type(Areachange),ALLOCATABLE::Arch_cl(:)    !areachange class 设定渐变段类型
    
    !/////倒虹吸类型////////////////
    Type::Invertedsiphon
        INTEGER::Serial_NUM !局部编号
        CHARACTER(len=12)::Name !节点名称
        REAL::Incoord !入口里程
        REAL::Outcoord !出口里程
        INTEGER::ParallelNum    !并排个数
        REAL::BottomWidth   !宽度
        REAL::Height    !高度
        REAL::Length    !长度
        REAL::InInvertElev  !入口底高程
        REAL::OutInvertElev  !出口底高程
        REAL::Roughness !倒虹吸糙率
        REAL::Insi_coeff !局部水力损失系数
    End type Invertedsiphon
    !//////////////////////////////////////////////
    Type(Invertedsiphon),ALLOCATABLE::Insi_cl(:)    !Invertedsiphon class 设定倒虹吸类型
    
    !////桥墩类型//////////
    Type::Bridge
        INTEGER::Serial_NUM !局部编号
        CHARACTER(len=12)::Name !节点名称
        REAL::Incoord !入口里程
        REAL::Outcoord !出口里程
        REAL::BottomWidth   !底宽
        REAL::SideSlopes    !边坡系数
        REAL::InInvertElev  !入口底高程
        REAL::OutInvertElev  !出口底高程
        REAL::Roughness !渠道糙率
        REAL::HinderWidth !总阻水宽度
    End type Bridge
    !////////////////////////////////////////
    Type(Bridge),ALLOCATABLE::Brid_cl(:)   !Bridge class设定渠道类型
    
    !/////节制闸类型////////////////
    Type::Controlgate
        INTEGER::Serial_NUM !局部编号
        CHARACTER(len=12)::Name !节点名称
        REAL::Incoord !入口里程
        REAL::Outcoord !出口里程
        REAL::Initial_gate_Zup  !初始水位
        INTEGER::ParallelNum    !并排个数
        REAL::SingleWidth   !宽度
        REAL::InvertElev  !底高程
        REAL::Opendegree    !闸门开度
        INTEGER:: RUNNum    !开启孔数
        REAL::G_A,G_B,G_C      !过闸流量系数与闸门开度之间的关系曲线
        REAL::Discoeff      !闸门过流系数
        REAL::Lenth
    End type Controlgate
    !//////////////////////////////////////////////
    Type(Controlgate),ALLOCATABLE::Cogate_cl(:)    !Controlgate class 设定节制闸类型
!    REAL,ALLOCATABLE::Opendegree(:,:)   !闸门开度过程
!    REAL,ALLOCATABLE::C_Opendegree(:,:)   !计算过程使用闸门开度

    !/////分水口类型////////////////
    Type::SideOutFlow
        INTEGER::Serial_NUM !局部编号
        CHARACTER(len=12)::Name !节点名称
        REAL::Incoord !入口里程
        REAL::Outcoord !出口里程
        REAL::OutDischarge  !分水流量
    End type SideOutFlow
    !//////////////////////////////////////////////
    Type(SideOutFlow),ALLOCATABLE::SideOF_cl(:)    !SideOutFlow class 设定节制闸类型
    REAL,ALLOCATABLE::OutDischarge(:,:)     !分水口出流过程
    REAL,ALLOCATABLE::C_OutDischarge(:,:)     !计算过程使用分水口出流
    
    !/////泵站类型////////////////
    Type::Pumpingstation
        INTEGER::Serial_NUM !局部编号
        CHARACTER(len=12)::Name !节点名称
        REAL::Incoord !入口里程
        REAL::Outcoord !出口里程
        REAL::Initial_pump_Zup  !初始水位
        REAL::InvertElev  !底高程
        INTEGER::UNIT_NUM   !泵站机组个数
        INTEGER,ALLOCATABLE::Pump_Type(:)   !泵站中水泵类型的编号
        REAL,ALLOCATABLE::Coefficient_A(:)   !二次项系数序列
        REAL,ALLOCATABLE::Coefficient_B(:)    !一次项系数序列
        REAL,ALLOCATABLE::Coefficient_C(:)    !常数项序列
    End type Pumpingstation
    !//////////////////////////////////////////////
    Type(Pumpingstation),ALLOCATABLE::PumpST_cl(:)    !Pumpingstation class 设定泵站类型
    Type(Pumpingstation),ALLOCATABLE::Ini_PumpST_cl(:)
    
    !/////水泵类型////////////////
    Type::Pump
        INTEGER::Serial_NUM !类型编号
        CHARACTER(len=12)::Name !水泵名称
        INTEGER::Coefficient_NUM   !转角个数
        REAL,ALLOCATABLE::Bladeangle_Serial(:)   !转角序列
        REAL,ALLOCATABLE::Coefficient_A(:)   !二次项系数序列
        REAL,ALLOCATABLE::Coefficient_B(:)    !一次项系数序列
        REAL,ALLOCATABLE::Coefficient_C(:)    !常数项序列
    END TYPE Pump
    !//////////////////////////////////////////////
    Type(Pump),ALLOCATABLE::Pump_cl(:)    !Pumpingstation class 设定水泵类型
    
    !/////转角过程时间序列//////////
    Type::Bladeangle_TimeSerial
        REAL,ALLOCATABLE::Each_PS(:,:)   !泵站中每个泵的转角过程   （转角，时间）
    END TYPE Bladeangle_TimeSerial
    Type(Bladeangle_TimeSerial),ALLOCATABLE::Bladeangle_cl(:)       !括号里是泵站编号
    Type(Bladeangle_TimeSerial),ALLOCATABLE::C_Bladeangle_cl(:)       !括号里是泵站编号
    
    !//////最大循环次数
    INTEGER::Cycle_NUM  !每一迭代步的最大循环数
    !///////线性方程组中的未知量
    
    !/////收敛条件
    INTEGER::endflag=20  !5代连续满足条件则结束
    
    !/////开始条件
    INTEGER::startflag=1
    
    !/////最佳输出条件
    INTEGER::bestoutputflag=1
    
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


    !///////方程个数
    INTEGER::Eq_NUM 
    !///////定义系数矩阵

    !////同化观测量
    
    !///////////数据同化集合滤波的相关变量说明
    INTEGER::NUMMember !集合中成员个数

    INTEGER::Boundary_sign_UP,Boundary_sign_DOWN  !边界条件控制符，0-流量，1-水位，2-水位~流量关系
    INTEGER::Output_sign    !输出结果控制符，0-只输出节制闸处结果；1-输出各元件断面结果
    REAL::m_coeff    !末端闸门过闸流量系数计算
    INTEGER::Control_sign  !同化模型和模拟模型的控制符，模拟为0，同化为1
    INTEGER::Statetype   !状态量类型（0-糙率；1-局损；2-糙率+局损）
    INTEGER::Obsertype   !观测量类型（0-水位；1-流量；2-水位+流量）
    REAL::PPQ, VVQ  !水位~流量关系断面系数
    


    !////////优化过程使用变量
    !////////泵站
    REAL::Coefficient_1_5(34,4) !存储第1-5级泵站的转角系数,已经给定，第一维给编号，第二-四维给系数，分别对应A,B,C；
    INTEGER::Bladeangle_1_5(34,3) !每一组转角系数对应的转角，第一-三维分别表示1、2、3号机组的转角值
    REAL::Coefficient_6(55,4)   !存储第6级泵站的转角系数，已经给定，第一维给编号，第二-四维给系数，分别对应A,B,C；
    INTEGER::Bladeangle_6(55,3)     !每一组转角系数对应的转角，第一-三维分别表示1、2、3号机组的转角值
    REAL,ALLOCATABLE::Forebay_MAX(:),Forebay_MIN(:),Outlet_MAX(:),Outlet_MIN(:),DESIGN_ZIN(:) !存储各级泵站站前、站后最高、最低水位约束|   
    INTEGER::Bladeangle_varibale(30)

    !节制闸
!    REAL,ALLOCATABLE::ST_STATE(:,:)     !存储各节制闸每小时闸前水位、闸后水位
!    REAL,ALLOCATABLE::DS_STATE(:)
!    INTEGER,ALLOCATABLE::SG_POSTION(:)	!节制闸位置
!     INTEGER,ALLOCATABLE::DM_POSTION(:)  !节制闸+首尾断面位置
!     INTEGER,ALLOCATABLE::ELE_POSTION(:) !渠道、桥墩以外的元件位置
!    REAL,ALLOCATABLE::GAV_Opendegree(:)
!    REAL,ALLOCATABLE::Deltadegree(:,:)
    REAL,ALLOCATABLE::ZQSW_TARGET(:)    !节制闸闸前目标水位，有几个节制闸就有几个
    REAL::D_WL_TARGET   !最下游边界处的目标水位
    REAL,ALLOCATABLE::ST_CON_UP(:)  !水位约束上限
    REAL,ALLOCATABLE::ST_CON_DOWN(:)    !水位约束下限
    REAL::DS_CON_UP,DS_CON_DOWN !下游水位边界，上限、下限
    REAL,ALLOCATABLE::h_new_INIT(:),QQ_INIT(:,:)    !初始水位、初始流量
    REAL,ALLOCATABLE::Opendegree_INIT(:,:)  !初始开度
    
    
    REAL,ALLOCATABLE::BestOpendegree(:,:)
    REAL,ALLOCATABLE::BestST_STATE(:,:)
    REAL,ALLOCATABLE::BestQControlGate(:,:)
    REAL,ALLOCATABLE::BestVolume(:,:)
End


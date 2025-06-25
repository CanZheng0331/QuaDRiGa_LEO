clc;
clear;
close all;


%% 卫星参数设置
satant = qd_arrayant('parabolic', 0.25, 20e9, [], 5, 1, [], 38.5); % LEO卫星天线（KA频段下行）
ueant  = qd_arrayant('parabolic', 0.6,  20e9, [], 5, 1, [], 39.7); % 用户终端天线（VSAT）

%% 卫星轨道位置
leo_sat = qd_satellite('custom', qd_satellite.R_e + 600, 0, 63.4, -28.8, 44.55, 0);
leo_sat = leo_sat.init_tracks([-5, 38.85]);  % 参考位置：西班牙中部
leo_sat.calc_orientation([], -90*pi/180);    % 天线指向正下方
leo_sat.name = 'LEO600';

%% 创建仿真布局
l = qd_layout;                                % 初始化布局
l.simpar.center_frequency = 20e9;             % 中心频率20GHz
l.simpar.show_progress_bars = 0;              % 关闭进度条

% 添加一个卫星（一个波束）
l.no_tx = 1;                                  % 单个发射机（卫星）
l.tx_track(1,1) = copy(leo_sat);              % 复制卫星轨道
l.tx_track(1,1).name = 'LEOBeam';            % 波束名称
l.tx_array(1,1) = satant;                     % 卫星天线

% 添加10个用户
l.no_rx = 10;                                 % 10个用户终端
beam_radius = 8e3;                            % 波束覆盖半径8km
l.randomize_rx_positions(beam_radius, 0, 0, 0); % 在16km×16km区域内随机分布用户

% 设置用户天线指向卫星
tx_pos = l.tx_track(1,1).initial_position;    % 卫星位置
for ir = 1:l.no_rx
    rx_pos = l.rx_track(1,ir).initial_position; % 用户位置
    rt = tx_pos - rx_pos;                      % 指向向量
    rt = rt/norm(rt);                          % 单位化
    
    % 正确构建3×1的朝向向量 [横滚角; 俯仰角; 偏航角]
    orientation = zeros(3, 1);                 % 初始化朝向向量
    orientation(2, 1) = asin(rt(3));           % 俯仰角计算 (Y轴旋转)
    orientation(3, 1) = atan2(rt(2), rt(1));   % 偏航角计算 (Z轴旋转)
    
    l.rx_track(1,ir).orientation = orientation; % 设置天线方向
    l.rx_array(1,ir) = ueant;                  % 用户天线
end

%% 设置传播场景
l.set_scenario('QuaDRiGa_NTN_Rural_LOS');     % NTN农村LOS场景

%% 显示布局
%figure;
l.visualize([],[],0);                         % 可视化布局
hold on;

% 获取卫星的位置
tx_pos = l.tx_track(1,1).initial_position;    % 卫星位置

% 绘制波束覆盖范围（圆形区域）
theta = linspace(0, 2*pi, 100);               % 用于生成圆的角度
x = beam_radius * cos(theta) + tx_pos(1);     % 圆的X坐标
y = beam_radius * sin(theta) + tx_pos(2);     % 圆的Y坐标
plot(x, y, 'r--', 'LineWidth', 2);             % 用红色虚线绘制波束范围

% 添加标题和标签
title('LEO卫星与用户分布');
xlabel('东向距离 (m)');
ylabel('北向距离 (m)');
grid on;

%% 运行仿真
c = l.get_channels;                           % 获取信道系数

%% 计算性能指标
PTx = 21.5;                % 发射功率(dBm)
NF = 1.2;                  % 噪声系数(dB)
BW = 400e6;                % 带宽(400MHz)
noise_thermal = -228.6 + 10*log10(290) + 30; % 热噪声(dBm/Hz)
noise_power = noise_thermal + 10*log10(BW) + NF; % 总噪声功率(dBm)

CPL = zeros(1, l.no_rx);   % 初始化耦合损耗
SINR_dB = zeros(1, l.no_rx); % 初始化SINR

for ir = 1:l.no_rx
    % 计算接收功率|H|²
    ch_coeff = c(ir,1).coeff;              % 获取信道系数
    power = sum(abs(ch_coeff(:)).^2);       % 计算信道功率
    
    % 计算耦合损耗(dB)
    CPL(ir) = -10*log10(power);
    
    % 计算接收功率(dBm)
    Prx_dBm = PTx + 10*log10(power);
    
    % 计算SINR(dB)
    SINR_dB(ir) = Prx_dBm - noise_power;
end

%% 显示结果
disp('=== 性能指标汇总 ===');
disp(['卫星类型: ', leo_sat.name]);
disp(['中心频率: ', num2str(l.simpar.center_frequency/1e9), ' GHz']);
disp(['发射功率: ', num2str(PTx), ' dBm']);
disp(['噪声系数: ', num2str(NF), ' dB']);
disp(['接收带宽: ', num2str(BW/1e6), ' MHz']);
fprintf('\n');

disp('用户耦合损耗 (dB):');
disp(CPL);
disp(['平均耦合损耗: ', num2str(mean(CPL)), ' dB']);

fprintf('\n');
disp('用户SINR (dB):');
disp(SINR_dB);
disp(['平均SINR: ', num2str(mean(SINR_dB)), ' dB']);

%% 绘制结果
figure;
plot(1:10, CPL, 'o-');
title('用户耦合损耗');
xlabel('用户编号'); ylabel('损耗 (dB)');
grid on;

figure;
plot(1:10, SINR_dB, 's-');
title('用户信噪比 (SINR)');
xlabel('用户编号'); ylabel('SINR (dB)');
grid on;
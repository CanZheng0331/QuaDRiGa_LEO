clc;
clear;
close all;

%% 卫星参数设置
% 卫星天线：4x4均匀平面阵列（UPA），中心频率1.9925 GHz
satant = qd_arrayant('custom', 4, 4, 1.9925e9, [], [], [], 38.5); % 4x4 UPA，增益38.5 dBi
ueant = qd_arrayant('omni'); % 用户终端：全向天线，适合手持终端

%% 卫星轨道位置
leo_sat = qd_satellite('custom', qd_satellite.R_e + 350e3, 0, 53, -28.8, 44.55, 0); % 轨道高度350 km，倾角53°
leo_sat = leo_sat.init_tracks([-5, 38.85]); % 参考位置：西班牙中部
leo_sat.calc_orientation([], -90*pi/180); % 天线指向正下方
leo_sat.name = 'StarlinkD2C';

%% 创建仿真布局
l = qd_layout; % 初始化布局
l.simpar.center_frequency = 1.9925e9;   % 中心频率1.9925 GHz
l.simpar.show_progress_bars = 0;        % 关闭进度条

% 添加一个卫星（一个波束）
l.no_tx = 1; % 单个发射机（卫星）
l.tx_track(1,1) = copy(leo_sat); % 复制卫星轨道
l.tx_track(1,1).name = 'D2CBeam'; % 波束名称
l.tx_array(1,1) = satant; % 卫星阵列天线

% 添加10个用户
l.no_rx = 50; % 10个用户
beam_radius = 50e3; % 波束覆盖半径50 km

% 获取卫星位置和地面投影点
tx_pos = l.tx_track(1).initial_position; % 卫星位置 (x, y, z)
ground_pos = [tx_pos(1); tx_pos(2); 0]; % 地面投影点 (x, y, z=0)，3×1列向量

% 随机生成用户位置并平移到地面投影点
l.randomize_rx_positions(beam_radius, 0, 0, 0); % 生成随机位置（100km×100km）
for ir = 1:l.no_rx
    rx_pos = l.rx_track(1,ir).initial_position; % 获取用户原始位置 (3×1)
    rx_pos = [rx_pos(1) + ground_pos(1); rx_pos(2) + ground_pos(2); 0]; % 平移，保持3×1
    l.rx_track(1,ir).initial_position = rx_pos; % 更新用户位置
    
    rt = tx_pos - rx_pos; % 指向向量
    rt = rt/norm(rt); % 单位化
    
    % 构建3×1朝向向量 [横滚角; 俯仰角; 偏航角]
    orientation = zeros(3, 1);
    orientation(2,1) = asin(rt(3)); % 俯仰角 (Y轴旋转)
    orientation(3,1) = atan2(rt(2), rt(1)); % 偏航角 (Z轴旋转)
    
    l.rx_track(1,ir).orientation = orientation; % 设置方向
    l.rx_array(1,ir) = ueant; % 用户全向天线
end

% 调试输出：验证坐标和距离
disp('Satellite position (x, y, z):');
disp(tx_pos');
disp('Ground projection (x, y, z):');
disp(ground_pos');
disp('User positions (x, y, z):');
for ir = 1:l.no_rx
    rx_pos = l.rx_track(1,ir).initial_position;
    disp(rx_pos');
    % 计算用户到地面投影点的距离
    dist = norm(rx_pos(1:2) - ground_pos(1:2));
    disp(['Distance from user ', num2str(ir), ' to ground_pos: ', num2str(dist), ' m']);
end

%% 设置传播场景
l.set_scenario('MIMOSA_10-45_LOS'); % MIMOSA场景，适合S波段卫星通信

%% 显示布局
l.visualize([],[],0); % 可视化布局
hold on;

% 绘制地面波束覆盖范围
theta = linspace(0, 2*pi, 100); % 生成圆的角度
x = beam_radius * cos(theta) + ground_pos(1); % 圆的X坐标
y = beam_radius * sin(theta) + ground_pos(2); % 圆的Y坐标
plot(x, y, 'r--', 'LineWidth', 2); % 红色虚线波束范围

% 添加用户编号标签
for ir = 1:l.no_rx
    rx_pos = l.rx_track(1,ir).initial_position;
    text(rx_pos(1) + 1000, rx_pos(2) + 1000, num2str(ir), ...
        'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
end

% 调整坐标刻度，聚焦波束覆盖区域
xlim([ground_pos(1) - beam_radius*1.2, ground_pos(1) + beam_radius*1.2]); % x轴范围±60 km
ylim([ground_pos(2) - beam_radius*1.2, ground_pos(2) + beam_radius*1.2]); % y轴范围±60 km

% 添加标题和标签
title('Starlink D2C卫星与用户分布');
xlabel('东向距离 (m)');
ylabel('北向距离 (m)');
grid on;
axis equal; % 确保x、y轴比例相同

%% 运行仿真
c = l.get_channels; % 获取信道系数

%% 计算性能指标
PTx = 21.5; % 发射功率(dBm)
NF = 1.2; % 噪声系数(dB)
BW = 400e6; % 带宽(400MHz)
noise_thermal = -228.6 + 10*log10(290) + 30; % 热噪声(dBm/Hz)
noise_power = noise_thermal + 10*log10(BW) + NF; % 总噪声功率(dBm)

CPL = zeros(1, l.no_rx); % 初始化耦合损耗
SINR_dB = zeros(1, l.no_rx); % 初始化SINR

for ir = 1:l.no_rx
    % 计算接收功率|H|²
    ch_coeff = c(ir,1).coeff; % 获取信道系数
    power = sum(abs(ch_coeff(:)).^2); % 计算信道功率
    
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
plot(1:50, CPL, 'o-');
title('用户耦合损耗');
xlabel('用户编号'); ylabel('损耗 (dB)');
grid on;

figure;
plot(1:50, SINR_dB, 's-');
title('用户信噪比 (SINR)');
xlabel('用户编号'); ylabel('SINR (dB)');
grid on;
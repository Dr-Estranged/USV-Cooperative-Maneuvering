%% run_simulation_4_cases.m
% 基准控制器测试：验证 4 种不同扰动下的系统表现
clear; clc; close all;

%% 1. 初始化参数
rng(42); % 固定随机数种子，保证"随机扰动"每次跑出来的图一模一样
M_matrix = diag([25.8, 33.8, 2.76]); 
model_params.M = M_matrix;           
model_params.M_inv = inv(M_matrix);  

ctrl_params.K1 = diag([0.5, 0.5, 0.5]);    
ctrl_params.K2 = diag([5.0, 5.0, 5.0]);    
ctrl_params.zeta1 = [4; 4; 4];    
ctrl_params.zeta2 = [1.5; 1.5; 1.5]; 

dt = 0.01; 
T_end = 50;
time = 0:dt:T_end;
N = length(time);

% 预定义 4 种工况名称，方便画图调用
case_names = {'(a) 常数型扰动', '(b) 周期型扰动', '(c) 随机扰动', '(d) 状态依赖型扰动'};
results = cell(4, 1); % 用一个元胞数组把4次实验的数据全存下来

disp('开始进行 4 种扰动工况的批量仿真...');

%% 2. 核心大循环 (1到4)
for case_id = 1:4
    fprintf('正在运行工况 %d/4: %s...\n', case_id, case_names{case_id});
    
    % [核心切换] 告诉物理引擎当前处于哪种扰动模式
    model_params.dist_type = case_id;
    
    % 重置状态
    X = zeros(6, 1);       
    X(1:2) = [6; 1];       
    X(3) = pi/2;           
    filter_state = zeros(6, 1); 
    prev_psi_d = 0;
    
    % 临时数据记录矩阵
    log_X = zeros(N, 6);
    log_z1 = zeros(N, 3);
    log_tau = zeros(N, 3);
    
    % 跑单次仿真
    for k = 1:N
        t = time(k);
        
        % 轨迹生成 (8字形)
        w1 = 0.1; w2 = 0.2;
        x_d = (10 + 2*sin(0.5*w1*t)) * cos(w1*t);
        y_d = (8 + 3*cos(0.4*w1*t)) * sin(w2*t);
        
        u_x = 10 + 2*sin(0.5*w1*t); u_x_dot = 1.0*w1*cos(0.5*w1*t);
        v_x = cos(w1*t);            v_x_dot = -w1*sin(w1*t);
        dx_d = u_x_dot * v_x + u_x * v_x_dot;
        
        u_y = 8 + 3*cos(0.4*w1*t);  u_y_dot = -1.2*w1*sin(0.4*w1*t);
        v_y = sin(w2*t);            v_y_dot = w2*cos(w2*t);
        dy_d = u_y_dot * v_y + u_y * v_y_dot;
        
        psi_d = atan2(dy_d, dx_d);
        if k == 1
            r_d = 0; prev_psi_d = psi_d;
        else
            d_psi = psi_d - prev_psi_d;
            while d_psi > pi, d_psi = d_psi - 2*pi; end
            while d_psi < -pi, d_psi = d_psi + 2*pi; end
            r_d = d_psi / dt; prev_psi_d = psi_d;
        end
        
        eta_bar = [x_d; y_d; psi_d];
        eta_dot_bar = [dx_d; dy_d; r_d]; 
        
        % 控制与环境更新
        [tau, debug, next_filter_state] = isv_controller_exact(X, eta_bar, eta_dot_bar, filter_state, model_params, ctrl_params, dt);
        dX = isv_plant(t, X, tau, model_params);
        
        X = X + dX * dt;
        filter_state = next_filter_state;
        
        % 记录
        log_X(k,:) = X';
        log_z1(k,:) = debug.z1';
        log_tau(k,:) = tau';
    end
    
    % 跑完一次，把所有记录打包存入 results 元胞中
    results{case_id}.log_X = log_X;
    results{case_id}.log_z1 = log_z1;
    results{case_id}.log_tau = log_tau;
end
disp('批量仿真结束！正在生成高规格对比图...');

%% 3. 画图排版：2x2 子图矩阵法
set(0, 'DefaultTextInterpreter', 'tex');
color_isv = [94, 181, 229]/255;
color_ref = [169, 84, 156]/255;

% 提前算出画轨迹用的参考路径
x_ref_plot = (10 + 2*sin(0.5*w1*time)) .* cos(w1*time);
y_ref_plot = (8 + 3*cos(0.4*w1*time)) .* sin(w2*time);

% === 大图 1: 四种扰动下的轨迹跟踪对比 ===
fig1 = figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [50, 50, 1000, 800]);

for i = 1:4
    subplot(2, 2, i);
    plot(x_ref_plot, y_ref_plot, '--', 'Color', color_ref, 'LineWidth', 1.5); hold on;
    plot(results{i}.log_X(:,1), results{i}.log_X(:,2), '-', 'Color', color_isv, 'LineWidth', 1.5);
    
    plot(results{i}.log_X(1,1), results{i}.log_X(1,2), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', color_isv, 'LineWidth', 1.5);
    plot(results{i}.log_X(end,1), results{i}.log_X(end,2), 's', 'MarkerSize', 6, 'MarkerFaceColor', color_isv);
    
    grid on; axis equal; box on;
    title(['\textbf{', case_names{i}, '}'], 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('x (m)'); ylabel('y (m)');
    if i == 1 % 只在第一个子图画图例
        legend({'参考路径', '实际轨迹'}, 'Location', 'best');
    end
end
sgtitle('\textbf{四种扰动下的轨迹跟踪表现 (基准控制器)}', 'Interpreter', 'latex', 'FontSize', 16);

% === 大图 2: 四种扰动下的跟踪误差对比 ===
fig2 = figure(2); clf;
set(gcf, 'Color', 'w', 'Position', [1100, 50, 1000, 800]);

for i = 1:4
    subplot(2, 2, i);
    plot(time, results{i}.log_z1(:,1), 'LineWidth', 1.5, 'Color', [33,158,188]/255); hold on;
    plot(time, results{i}.log_z1(:,2), 'LineWidth', 1.5, 'Color', [250,134,000]/255);
    plot(time, results{i}.log_z1(:,3), 'LineWidth', 1.5, 'Color', [2,48,74]/255);
    
    grid on; xlim([0, T_end]); 
    ylim([-2.5, 4]); % 统一 Y 轴量程，方便对比！
    
    title(['\textbf{', case_names{i}, '}'], 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('时间 (s)'); ylabel('误差幅值');
    if i == 1
        legend({'纵向误差 x_e', '横向误差 y_e', '航向误差 \psi_e'}, 'Location', 'northeast');
    end
end
sgtitle('\textbf{四种扰动下的船体坐标系跟踪误差对比}', 'Interpreter', 'latex', 'FontSize', 16);

%% 4. 导出高清 PDF
set(findall(fig1, '-property', 'FontName'), 'FontName', 'Microsoft YaHei');
set(findall(fig2, '-property', 'FontName'), 'FontName', 'Microsoft YaHei');
try
    exportgraphics(fig1, '4_Cases_Trajectory.pdf');
    exportgraphics(fig2, '4_Cases_Errors.pdf');
    disp('PDF 导出成功！');
catch
    disp('自动导出 PDF 失败，请手动保存。');
end
%% run_multi_agent_simulation.m
% 复现论文: Noncooperative Game-based Cooperative Maneuvering of ISVs
% 重点复现 Action Estimators 部分 (Section III-A) 和 Fig 3-7
clear; clc; close all;

%% 1. 仿真参数设置 (对应论文 Section V)
M = 5;              % ISV数量
dt = 0.01;          % 时间步长
T_end = 400;        % 仿真总时长
time = 0:dt:T_end;
N_steps = length(time);

% 通信拓扑邻接矩阵 (对应论文 Fig. 2)
A_matrix = [0 1 0 0 0; 
            1 0 1 0 1; 
            0 1 0 1 0; 
            0 0 1 0 1; 
            0 1 0 1 0];

% 参考路径参数 (论文公式中的 R_i)
R_list = [150, 170, 190, 210, 230];

%% 2. 物理模型参数 (对应论文 Section V)
M_matrix = diag([25.8, 33.8, 2.76]); 
model_params.M = M_matrix;
model_params.M_inv = inv(M_matrix);

%% 3. 控制器参数 (对应论文 Section V)
% 估计器参数 (公式 9, 10)
est_params.mu    = diag([5, 5, 5]);   
est_params.kappa = 5;                  

% 控制器参数 (公式 23, 24)
ctrl_params.K1 = diag([1, 1, 1]);     
ctrl_params.K2 = diag([10, 10, 10]);  
ctrl_params.zeta1 = [40; 40; 40];
ctrl_params.zeta2 = [0.05; 0.05; 0.05];

% Nash均衡寻优参数 (公式 34-37)
delta_eta   = diag([0.1, 0.1, 0.1]);  
delta_theta = 0.1;                     
k_p3        = 1;                       

%% 4. 状态初始化 - 关键：制造初始偏差以展示收敛过程
% 初始路径参数
theta_init_val = 0;

% ========== 关键修改：ISV初始位置偏离Nash均衡目标 ==========
X_all = zeros(6, M);
initial_offset_x = [30, -20, 40, -30, 25];   % x方向偏移
initial_offset_y = [20, 35, -25, 40, -30];   % y方向偏移

for i = 1:M 
    % Nash均衡目标位置 (参考路径上的点)
    nash_x = R_list(i) * cos(theta_init_val);
    nash_y = R_list(i) * sin(theta_init_val);
    nash_psi = 1.03 * theta_init_val + 0.4855 * pi;
    
    % ISV实际初始位置 = Nash目标 + 偏移 (制造偏差)
    X_all(1, i) = nash_x + initial_offset_x(i);  
    X_all(2, i) = nash_y + initial_offset_y(i);  
    X_all(3, i) = nash_psi + 0.1*(rand()-0.5);   
    X_all(4:6, i) = 0;  
end 

% 路径参数 θ 初始化 (各ISV有不同的初始θ，制造不一致)
Theta = [0.05, 0.02, 0.08, 0.01, 0.06];  % 不同初始值

% Nash均衡目标初始化 (从参考路径上的点开始)
Nash_Eta = zeros(3, M);
for i = 1:M
    Nash_Eta(1, i) = R_list(i) * cos(Theta(i)); 
    Nash_Eta(2, i) = R_list(i) * sin(Theta(i));
    Nash_Eta(3, i) = 1.03 * Theta(i) + 0.4855 * pi;
end
Nash_Theta = Theta; 

% 估计器初始化 (初始时有估计误差)
Eta_bar_hat = zeros(3, M, M); 
Theta_bar_hat = zeros(1, M, M);

for i = 1:M
    for j = 1:M
        if i == j
            Eta_bar_hat(:, j, i) = Nash_Eta(:, j); 
            Theta_bar_hat(1, j, i) = Nash_Theta(j);
        else
            % 对其他ISV的估计有初始误差
            Eta_bar_hat(:, j, i) = Nash_Eta(:, j) + 10*(rand(3,1)-0.5);
            Theta_bar_hat(1, j, i) = Nash_Theta(j) + 0.05*(rand()-0.5);
        end
    end
end

% 滤波器状态初始化
filter_states = zeros(6, M); 

%% 5. 数据记录数组
log_X = zeros(N_steps, 6, M);           
log_Theta = zeros(N_steps, M);          
log_Nash_Eta = zeros(N_steps, 3, M);    
log_Nash_Theta = zeros(N_steps, M);     
log_Est_Error_Eta = zeros(N_steps, M);  
log_Est_Error_Theta = zeros(N_steps, M);
log_Track_Error = zeros(N_steps, 3, M); 

%% 6. 主仿真循环
fprintf('开始仿真... 总时长: %.1f s\n', T_end);
tic;

for k = 1:N_steps
    t = time(k);
    
    %% Step 1: 计算干扰 (用于完美补偿)
    Disturbances = zeros(3, M);
    for i = 1:M
        nu = X_all(4:6, i);
        D_lin = diag([50, 50, 10]); 
        D_nonlin = diag([10*abs(nu(1)), 10*abs(nu(2)), 5*abs(nu(3))]);
        f_nu = -(D_lin + D_nonlin) * nu; 
        w_t = [10 * sin(0.1*t); 10 * cos(0.1*t); 1 * sin(0.05*t)];
        Disturbances(:, i) = f_nu + w_t;
    end
    
    %% Step 2: Action Estimators 更新 (公式 9, 10)
    [d_Eta_bar_hat, d_Theta_bar_hat] = update_estimators(...
        Eta_bar_hat, Theta_bar_hat, Nash_Eta, Nash_Theta, ...
        A_matrix, est_params.mu, est_params.kappa);
    
    %% Step 3: Nash均衡寻优 (公式 34-37)
    d_Nash_Eta   = zeros(3, M);
    d_Nash_Theta = zeros(1, M);
    d_Theta      = zeros(1, M);
    
    for i = 1:M
        My_Est_Eta_bar = squeeze(Eta_bar_hat(:, :, i)); 
        My_Est_Theta_bar = squeeze(Theta_bar_hat(:, :, i))'; 
        
        my_eta_bar = Nash_Eta(:, i);
        my_theta_bar = Nash_Theta(i);
        
        [grad_J1, grad_J2] = calculate_gradients(i, my_eta_bar, my_theta_bar, ...
                                                  My_Est_Eta_bar, My_Est_Theta_bar, t, R_list);
        
        d_Nash_Eta(:, i) = -delta_eta * grad_J1;
        d_Nash_Theta(i)  = -delta_theta * grad_J2;
        
        z3 = Theta(i) - Nash_Theta(i);
        d_Theta(i) = -k_p3 * z3 + d_Nash_Theta(i);
    end
    
    %% Step 4: 控制器计算 (公式 23, 24)
    Tau_cmd = zeros(3, M);
    for i = 1:M
        target_eta = Nash_Eta(:, i);
        target_dot = d_Nash_Eta(:, i);
        [tau_i, ~, next_filt] = isv_controller_exact(X_all(:,i), target_eta, target_dot, ...
                                    filter_states(:,i), model_params, ctrl_params, dt, Disturbances(:,i));
        Tau_cmd(:, i) = tau_i;
        filter_states(:, i) = next_filt;
    end
    
    %% Step 5: 物理模型更新
    dX_all = zeros(6, M);
    for i = 1:M
        dX_all(:, i) = isv_plant(t, X_all(:,i), Tau_cmd(:,i), model_params);
    end
    
    %% Step 6: 欧拉积分
    X_all         = X_all + dX_all * dt;
    Eta_bar_hat   = Eta_bar_hat + d_Eta_bar_hat * dt;
    Theta_bar_hat = Theta_bar_hat + d_Theta_bar_hat * dt;
    Nash_Eta      = Nash_Eta + d_Nash_Eta * dt;
    Nash_Theta    = Nash_Theta + d_Nash_Theta * dt;
    Theta         = Theta + d_Theta * dt;
    
    %% Step 7: 数据记录
    log_X(k, :, :) = X_all;
    log_Theta(k, :) = Theta;
    log_Nash_Eta(k, :, :) = Nash_Eta;
    log_Nash_Theta(k, :) = Nash_Theta;
    
    for i = 1:M
        eta_err_sum = 0;
        theta_err_sum = 0;
        for j = 1:M
            if j ~= i
                eta_err_sum = eta_err_sum + norm(Eta_bar_hat(:,j,i) - Nash_Eta(:,j));
                theta_err_sum = theta_err_sum + abs(Theta_bar_hat(1,j,i) - Nash_Theta(j));
            end
        end
        log_Est_Error_Eta(k, i) = eta_err_sum / (M-1);
        log_Est_Error_Theta(k, i) = theta_err_sum / (M-1);
    end
    
    for i = 1:M
        log_Track_Error(k, :, i) = X_all(1:3, i) - Nash_Eta(:, i);
    end
    
    if mod(k, round(N_steps/10)) == 0
        fprintf('进度: %.0f%%\n', k/N_steps*100);
    end
end

elapsed = toc;
fprintf('仿真完成! 耗时: %.2f s\n', elapsed);

%% 7. 绘图 (已根据论文排版修正放大框位置)
% 自定义颜色
colors = [216, 56, 58; 84, 179, 69; 248,172,140; 199, 109, 162; 40, 120, 181] / 255;

% 创建输出文件夹
script_path = fileparts(mfilename('fullpath'));
if isempty(script_path), script_path = pwd; end
output_folder = fullfile(script_path, 'figures');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
fprintf('图片将保存到: %s\n', output_folder);

%% 修改后的 Fig 3 绘图代码 (让它看起来像论文)
figure('Name', 'Fig 3: Reproduced Style', 'Position', [100, 100, 600, 500]);
hold on; axis equal; grid on;

% 2. 画 Nash 均衡轨迹 (论文中的 dotted line)
for i = 1:M
    % 画 Nash 均衡解 (你的代码里存了 log_Nash_Eta)
    plot(log_Nash_Eta(:, 1, i), log_Nash_Eta(:, 2, i), ':', 'Color', colors(i,:), 'LineWidth', 2);
end

% 3. 画 实际轨迹 (论文中的 dot dash / solid line)
h_traj = zeros(1, M);
for i = 1:M
    h_traj(i) = plot(log_X(:, 1, i), log_X(:, 2, i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % 画小船 (用简单的三角形或圆点代替)
    plot(log_X(end, 1, i), log_X(end, 2, i), 'd', ...
        'Color', 'k', 'MarkerFaceColor', colors(i,:), 'MarkerSize', 8);
end

xlabel('x (m)'); ylabel('y (m)');
title('Fig 3: Actual Trajectories vs Nash Equilibrium (Reproduced)');
legend(h_traj, {'ISV_1', 'ISV_2', 'ISV_3', 'ISV_4', 'ISV_5'}, 'Location', 'best');

exportgraphics(gcf, fullfile(output_folder, 'Fig3_traj.pdf'), 'ContentType', 'vector');

%% Fig 4: x坐标收敛 (修正：放大框移至左下角)
figure('Name', 'Fig 4: x-coordinate', 'Position', [720, 100, 650, 450]);
% === 主图 ===
axes('Position', [0.13 0.11 0.775 0.815]);
hold on; grid on;
h_lines = [];
legend_labels = {};
for i = 1:M
    h1 = plot(time, log_X(:, 1, i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    h2 = plot(time, log_Nash_Eta(:, 1, i), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    h_lines = [h_lines, h1, h2];
    legend_labels = [legend_labels, {sprintf('x_%d', i), sprintf('x_%d^*', i)}];
end
xlim([0, 400]); 
% Y轴范围根据数据自动调整，通常在 -250 到 300 之间
xlabel('Time (s)'); ylabel('x (m)');
title('Fig 4: Evolution of x position');
legend(h_lines, legend_labels, 'Location', 'northeast', 'NumColumns', 2, 'Interpreter', 'tex', 'FontSize', 8);

% === 局部放大图 (Inset) ===
% 【修改点】位置改到左下角 [left bottom width height]
% 已修改：从 [0.22 0.22 0.25 0.25] 改为 [0.16 0.16 0.25 0.25]
axes('Position', [0.16 0.16 0.25 0.25]); 
box on; hold on; grid on;
for i = 1:M
    plot(time, log_X(:, 1, i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    plot(time, log_Nash_Eta(:, 1, i), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
end
xlim([0, 12]);   
ylim([130, 270]); % 确保能看到所有曲线的起点
set(gca, 'FontSize', 8, 'Color', 'w');

% === 箭头 ===
% 【修改点】箭头从放大框顶部指向左上角的曲线起点
% 坐标为归一化坐标 [x_start, x_end], [y_start, y_end]
% 已微调箭头起点以匹配新的方框位置
annotation('arrow', [0.25 0.155], [0.42 0.72], 'Color', 'k', 'LineWidth', 1);

exportgraphics(gcf, fullfile(output_folder, 'Fig4_x_position.pdf'), 'ContentType', 'vector');

%% Fig 5: y坐标收敛 (修正：放大框移至下方居中)
figure('Name', 'Fig 5: y-coordinate', 'Position', [720, 550, 650, 450]);
% === 主图 ===
axes('Position', [0.13 0.11 0.775 0.815]);
hold on; grid on;
h_lines = [];
legend_labels = {};
for i = 1:M
    h1 = plot(time, log_X(:, 2, i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    h2 = plot(time, log_Nash_Eta(:, 2, i), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    h_lines = [h_lines, h1, h2];
    legend_labels = [legend_labels, {sprintf('y_%d', i), sprintf('y_%d^*', i)}];
end
xlim([0, 400]);
xlabel('Time (s)'); ylabel('y (m)');
title('Fig 5: Evolution of y position');
legend(h_lines, legend_labels, 'Location', 'northeast', 'NumColumns', 2, 'Interpreter', 'tex', 'FontSize', 8);

% === 局部放大图 (Inset) ===
% 【修改点】位置改到下方靠右一点，避开曲线拱起的下方
axes('Position', [0.38 0.2 0.25 0.25]); 
box on; hold on; grid on;
for i = 1:M
    plot(time, log_X(:, 2, i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    plot(time, log_Nash_Eta(:, 2, i), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
end
xlim([0, 12]);
ylim([-40, 60]); % 涵盖起始波动
set(gca, 'FontSize', 8, 'Color', 'w');

% === 箭头 ===
annotation('arrow', [0.38 0.16], [0.32 0.45], 'Color', 'k', 'LineWidth', 1);

exportgraphics(gcf, fullfile(output_folder, 'Fig5_y_position.pdf'), 'ContentType', 'vector');

%% Fig 6: 估计器误差 (保持不变)
figure('Name', 'Fig 6: Estimator Errors', 'Position', [100, 550, 600, 400]);
subplot(2,1,1);
hold on; grid on;
for i = 1:M
    plot(time, log_Est_Error_Eta(:, i), 'Color', colors(i,:), 'LineWidth', 1.5);
end
xlabel('Time (s)'); ylabel('Position Est. Error');
title('Fig 6(a): Position estimation errors');
xlim([0, 50]);
legend({'\eta_1', '\eta_2', '\eta_3', '\eta_4', '\eta_5'}, 'Location', 'best', 'Interpreter', 'tex');

subplot(2,1,2);
hold on; grid on;
for i = 1:M
    plot(time, log_Est_Error_Theta(:, i), 'Color', colors(i,:), 'LineWidth', 1.5);
end
xlabel('Time (s)'); ylabel('\theta Est. Error');
title('Fig 6(b): Parameter estimation errors');
xlim([0, 50]);
legend({'\theta_1', '\theta_2', '\theta_3', '\theta_4', '\theta_5'}, 'Location', 'best', 'Interpreter', 'tex');
exportgraphics(gcf, fullfile(output_folder, 'Fig6_Estimator_Errors.pdf'), 'ContentType', 'vector');

%% Fig 7: 路径参数 θ (修正：放大框在左上角)
figure('Name', 'Fig 7: Path Parameters', 'Position', [100, 300, 650, 450]);
% === 主图 ===
axes('Position', [0.13 0.11 0.775 0.815]);
hold on; grid on;
theta_ref_line = 0.0103 * time;
h_lines = [];
legend_labels = {};
for i = 1:M
    h1 = plot(time, log_Theta(:, i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    h_lines = [h_lines, h1];
    legend_labels = [legend_labels, {sprintf('\\theta_%d', i)}];
end
h_ref = plot(time, theta_ref_line, 'k--', 'LineWidth', 2);
h_lines = [h_lines, h_ref];
legend_labels = [legend_labels, {'\theta_r'}];
xlim([0, 400]);
xlabel('Time (s)'); ylabel('\theta (rad)');
title('Fig 7: Evolution of path parameters \theta_i');
legend(h_lines, legend_labels, 'Location', 'southeast', 'Interpreter', 'tex');

% === 局部放大图 (Inset) ===
% 位置在左上角
axes('Position', [0.2 0.55 0.25 0.25]); 
box on; hold on; grid on;
for i = 1:M
    plot(time, log_Theta(:, i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
end
plot(time, theta_ref_line, 'k--', 'LineWidth', 2);
xlim([0, 5]);
ylim([0, 0.15]);
set(gca, 'FontSize', 8, 'Color', 'w');

% === 箭头 ===
annotation('arrow', [0.23 0.15], [0.55 0.18], 'Color', 'k', 'LineWidth', 1);

exportgraphics(gcf, fullfile(output_folder, 'Fig7_Path_Parameters.pdf'), 'ContentType', 'vector');
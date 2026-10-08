%% 主仿真脚本
clear; clc; close all;

% --- 1. 定义模型参数 (Assumption) ---
% 质量矩阵 M (对角阵示例)
M = diag([25.8, 33.8, 2.76]); 
model_params.M_inv = inv(M); % 提前求逆

% --- 2. 定义初始条件 ---
% 初始位置 [0,0,0], 初始速度 [0,0,0]
X0 = zeros(6, 1); 

% --- 3. 定义控制输入 tau (测试用) ---
% 假设我们给一个恒定的推力和转矩 (开环测试)
% 对应截图中的绿色 tau_i
my_tau = [100; 0; 5]; % [纵向推力100N, 横向力0, 转向力矩5Nm]

% --- 4. 设置仿真时间 ---
tspan = [0 50]; % 跑50秒

% --- 5. 运行 ode45 ---
% 定义匿名函数，把 tau 和 params 传进去
ode_func = @(t, X) isv_plant(t, X, my_tau, model_params);
%匿名函数，其实要求解的就是X的六个变量，对于 ode45 来说，它看到的 ode_func 就变成了一个标准的二参数函数，它不需要知道内部包含了推力 tau 或质量 M，它只管传 t 和 X 进去，拿回导数。
[T, X_hist] = ode45(ode_func, tspan, X0);
%X_hist(k, :) 就代表在时刻 T(k) 时的 [x, y, \psi, u, v, r]。
% --- 5.5  回代计算导数项 ---
% ode45 只返回了积分后的状态 X，我们需要把 X 带回物理模型算出导数 dX
n_steps = length(T);dX_hist = zeros(n_steps, 6); % 初始化矩阵 [N x 6]
for i = 1:n_steps
    % 取出当前时刻的状态和时间
    t_now = T(i); X_now = X_hist(i, :)'; % 转置为列向量以符合函数输入要求
    % 调用物理模型函数，重新计算导数
    d_state = isv_plant(t_now, X_now, my_tau, model_params);
    dX_hist(i, :) = d_state';
end
% 现在 dX_hist 里面存的就是 [x_dot, y_dot, psi_dot, u_dot, v_dot, r_dot]

% --- 6. 绘图结果 ---
figure(1);
set(gcf, 'Color', 'w'); 
% === 图 1: 轨迹 ===
subplot(2,1,1);
deep_blue =[0.2039, 0.6902, 0.9412]; 
plot(X_hist(:,1), X_hist(:,2), ...
    'Color', deep_blue, ...       % 设置颜色
    'LineWidth', 2, ...           % 加粗线条
    'LineStyle', '-');            % 实线
% 添加起点标记 (用红色圆圈标出起点)
hold on;
plot(X_hist(1,1), X_hist(1,2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(X_hist(1,1)+0.5, X_hist(1,2), 'Start', 'FontSize', 10); % 标注起点
hold off;
title('\textbf{ISV Trajectory (x-y Plane)}', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
grid on;       % 打开网格
box on;        % 打开外边框
axis equal;    % 保证 xy 轴比例一致 (画圆必备)
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman'); % 设置刻度字体
grid on; box on;
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
color_u =[195, 105, 175]/255; % 用于 u
color_v = [242, 115, 158]/255; % 用于 v
color_r =[0.4510, 0.5882, 0.8627]; % 用于 r
color_x = [0.9647, 0.8000, 0.9059];      % 用于 x的颜色
color_y = [0.2745, 0.7843, 0.9412]; % 用于 y的颜色

%% === 图 2: 速度全家桶 (u, v, r 在一张图) ===
figure(2); set(gcf, 'Color', 'w');
title('\textbf{Velocity States Comparison}', 'Interpreter', 'latex', 'FontSize', 14);
grid on; hold on;
xlabel('Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 12);

% --- 左 Y 轴：画线速度 u 和 v (单位 m/s) ---
yyaxis left
p_u = plot(T, X_hist(:,4), 'Color', color_u, 'LineWidth', 2, 'LineStyle', '-'); 
hold on;
p_v = plot(T, X_hist(:,5), 'Color', color_v, 'LineWidth', 2, 'LineStyle', '-');
ylabel('Linear Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'ycolor', '[0.765,0.412, 0.686]'); % 设置左轴文字为黑色（或者选其中一个线条颜色）

% --- 右 Y 轴：画角速度 r (单位 rad/s) ---
yyaxis right
p_r = plot(T, X_hist(:,6), 'Color', color_r, 'LineWidth', 2, 'LineStyle', '--');
ylabel('Angular Velocity (rad/s)', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'ycolor',[0.4510, 0.5882, 0.8627]); % 设置右轴文字颜色与 r 一致，方便区分

legend([p_u, p_v, p_r], ...
       {'$u$ (surge)', '$v$ (sway)', '$r$ (yaw rate)'}, ...
       'Interpreter', 'latex', 'Location', 'best');
box on;

%% === 图 3: 位置随时间变化 (x-t 和 y-t 在一张图) ===
figure(3); set(gcf, 'Color', 'w');
plot(T, X_hist(:,1), 'Color', color_x, 'LineWidth', 2); hold on;
plot(T, X_hist(:,2), 'Color', color_y, 'LineWidth', 2);
title('\textbf{Position Evolution over Time}', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Position (m)', 'Interpreter', 'latex', 'FontSize', 12);
legend({'$x(t)$', '$y(t)$'}, 'Interpreter', 'latex', 'Location', 'best');
grid on; box on;
%% === 图 4: 三轴加速度对比 (合并显示) ===
figure(4); set(gcf, 'Color', 'w');
clf; grid on; hold on;
% 定义颜色
color_u =[0.2980, 0.4000, 0.6980];
color_v =[0.6784, 0.4157, 0.7176];
color_r =[0.9529, 0.4510, 0.6353];
% 绘三条曲线
% 注意：为了看清稳态时的震荡，我们这里画在一起
p1 = plot(T, dX_hist(:,4), 'Color', color_u, 'LineWidth', 1.5);
p2 = plot(T, dX_hist(:,5), 'Color', color_v, 'LineWidth', 1.5, 'LineStyle', '-');
p3 = plot(T, dX_hist(:,6), 'Color', color_r, 'LineWidth', 1.5, 'LineStyle', '-.');
% 美化
title('\textbf{Comparison of Accelerations}', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Magnitude ($m/s^2$ or $rad/s^2$)', 'Interpreter', 'latex', 'FontSize', 12);
% 添加图例
legend([p1, p2, p3], ...
       {'$\dot{u}$ (Surge Acc)', '$\dot{v}$ (Sway Acc)', '$\dot{r}$ (Yaw Acc)'}, ...
       'Interpreter', 'latex', 'Location', 'best');
% 提示：前几秒 u_dot 会有一个巨大的启动尖峰 (约 4 m/s^2)，
% 导致后面的震荡看起来很小。可以手动缩放一下 Y 轴看看细节：
ylim([0, 1]); % 限制 Y 轴范围，忽略启动时的尖峰，专注于稳态震荡
%% === 导出 PDF 文件 ===
exportgraphics(figure(1), '轨迹图.pdf', 'ContentType', 'vector');
exportgraphics(figure(2), 'uvr-t图像.pdf', 'ContentType', 'vector');
exportgraphics(figure(3), 'xy-t图像.pdf', 'ContentType', 'vector');
exportgraphics(figure(4), 'νi的导数（即加速度）图像.pdf', 'ContentType', 'vector');
disp('PDF 图片已成功导出！');
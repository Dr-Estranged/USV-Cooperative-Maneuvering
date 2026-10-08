function [tau, debug, next_filter_state] = isv_controller_exact(X, eta_bar, eta_dot_bar, filter_state, phys_params, ctrl_params, dt, known_disturbance)
    % ISV_CONTROLLER_EXACT 反步法控制器
    % 对应论文 Section III-C, 公式 (23) - (27)
    %
    % 功能: 控制ISV跟踪Nash均衡目标 eta_bar
    %
    % 输入:
    %   X: [6x1] 状态向量 [x; y; psi; u; v; r]
    %   eta_bar: [3x1] Nash均衡位置目标
    %   eta_dot_bar: [3x1] Nash均衡位置目标的导数
    %   filter_state: [6x1] 滤波器状态
    %   phys_params: 物理参数结构体
    %   ctrl_params: 控制器参数结构体
    %   dt: 时间步长
    %   known_disturbance: [3x1] 已知干扰 (用于补偿)
    %
    % 输出:
    %   tau: [3x1] 控制力/力矩
    %   debug: 调试信息结构体
    %   next_filter_state: [6x1] 下一时刻滤波器状态
    
    %% 0. 状态提取
    eta = X(1:3);  % 位置 [x; y; psi]
    nu  = X(4:6);  % 速度 [u; v; r]
    psi = X(3);    % 艏向角
    r   = X(6);    % 角速度 (psi_dot)
    
    %% 1. 构建旋转矩阵
    c = cos(psi); s = sin(psi);
    R = [c, -s, 0;
         s,  c, 0;
         0,  0, 1];
     
    % 旋转矩阵导数
    R_dot = [-s, -c, 0;
              c, -s, 0;
              0,  0, 0] * r; 
    
    %% 2. 运动学回路 (公式 23)
    % 位置误差
    eta_err = eta - eta_bar;
    % 角度归一化
    eta_err(3) = atan2(sin(eta_err(3)), cos(eta_err(3)));
    
    % 体坐标系误差
    z1 = R' * eta_err;
    
    % 虚拟控制律
    alpha = -ctrl_params.K1 * z1 + R' * eta_dot_bar;
    
    %% 3. 非线性指令滤波器 (公式 24)
    [alpha_f, alpha_dot_f, next_filter_state] = nonlinear_command_filter(alpha, filter_state, ctrl_params, dt);
    
    %% 4. 动力学回路
    z2 = nu - alpha_f;
    
    %% 5. 控制律计算 (公式 27)
    M = phys_params.M;
    
    % 控制力矩 (含干扰补偿)
    tau = M * (-ctrl_params.K2 * z2 + alpha_dot_f) - R' * z1 - known_disturbance;
    
    % 控制力矩限幅
    tau_max = [500; 500; 200];
    tau = max(min(tau, tau_max), -tau_max);
    
    %% 6. 调试信息
    debug.z1 = z1;
    debug.z2 = z2;
    debug.alpha = alpha;
    debug.alpha_f = alpha_f;
end

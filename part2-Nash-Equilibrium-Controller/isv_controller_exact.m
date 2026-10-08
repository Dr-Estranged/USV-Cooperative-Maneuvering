function [tau, debug, next_filter_state] = isv_controller_exact(X, eta_bar, eta_dot_bar, filter_state, phys_params, ctrl_params, dt)
    % Inputs:
    %   X: [x; y; psi; u; v; r]
    %   eta_bar: 参考轨迹 [x_d; y_d; psi_d]
    %   eta_dot_bar: 参考速度
    
    %% 0. 基础变量计算
    eta = X(1:3); % [x; y; psi]
    nu  = X(4:6); % [u; v; r]
    psi = X(3);
    r   = X(6);   % r = psi_dot
    
    % 旋转矩阵 R(psi)
    c = cos(psi); s = sin(psi);
    R = [c, -s, 0;
         s,  c, 0;
         0,  0, 1];
         
    % 旋转矩阵的导数 R_dot = dR/dpsi * psi_dot
    R_dot = [-s, -c, 0;
              c, -s, 0;
              0,  0, 0] * r; 
    
    %% --- Step 1: 运动学回路 (Exact Formula) ---
    
    % 1. 计算位置误差 z1 (Body Frame)
    % 原文: z_{i,1} = R^T * (eta - eta_bar)
    eta_err = eta - eta_bar;
    eta_err(3) = wrapToPi(eta_err(3)); % 角度归一化
    
    z1 = R' * eta_err;
    
    % 2. 计算虚拟控制律 alpha (Exact Formula)
    % 原文: alpha = -K1*z1 - R_dot^T * (eta - eta_bar) + R^T * eta_dot_bar
    alpha = -ctrl_params.K1 * z1 - R_dot' * eta_err + R' * eta_dot_bar;
    
    %% --- Step 1.5: 非线性指令滤波器 ---
    [alpha_f, alpha_dot_f, next_filter_state] = nonlinear_command_filter(alpha, filter_state, ctrl_params, dt);
    
    %% --- Step 2: 动力学回路 ---
    
    % 1. 计算速度误差 z2 = nu - alpha_f
    z2 = nu - alpha_f;
    
    %% --- Step 3: 力矩计算 (Tau) ---
    % 原文: u_i = M * (-K2*z2 + alpha_dot_f) - z1
    M = phys_params.M;
    tau = M * (-ctrl_params.K2 * z2 + alpha_dot_f) - z1;
    
    %% Debug Info
    debug.z1 = z1;
    debug.z2 = z2;
end


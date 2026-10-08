function [alpha_f, alpha_dot_f, next_state] = nonlinear_command_filter(alpha, current_state, params, dt)
    % NONLINEAR_COMMAND_FILTER 非线性指令滤波器 (Levant微分器)
    % 对应论文公式 (24), (25)
    %
    % 公式:
    %   dot_alpha_f = -zeta1 * sig(alpha_f - alpha)^0.5 + alpha_d_f
    %   dot_alpha_d_f = -zeta2 * sig(alpha_f - alpha)^0
    %
    % 其中 sig(x)^a = |x|^a * sign(x)
    %
    % 输入:
    %   alpha: [3x1] 虚拟控制输入
    %   current_state: [6x1] 滤波器状态 [alpha_f; alpha_d_f]
    %   params: 参数结构体
    %   dt: 时间步长
    %
    % 输出:
    %   alpha_f: [3x1] 滤波后的虚拟控制
    %   alpha_dot_f: [3x1] 虚拟控制的导数
    %   next_state: [6x1] 下一时刻滤波器状态
    
    % 状态解包
    x1 = current_state(1:3);  % alpha_f (滤波后的输出)
    x2 = current_state(4:6);  % alpha_d_f (导数估计)
    
    % 参数
    z1 = params.zeta1; 
    z2 = params.zeta2; 
    
    % 误差
    err = x1 - alpha;
    
    % 平滑符号函数 (避免高频震荡)
    phi = 0.1;  % 边界层厚度
    soft_sign = tanh(err ./ phi); 
    
    % sig函数计算
    % sig(x)^0.5 = |x|^0.5 * sign(x)
    % sig(x)^0 = sign(x)
    sig_05 = abs(err).^0.5 .* soft_sign;
    sig_0  = soft_sign;
    
    % 动力学方程
    dx1 = -z1 .* sig_05 + x2;
    dx2 = -z2 .* sig_0;
    
    % 欧拉积分
    x1_new = x1 + dx1 * dt;
    x2_new = x2 + dx2 * dt;
    
    % 输出
    alpha_f = x1_new;
    alpha_dot_f = dx1;  % dx1 就是 alpha_f 的导数
    
    next_state = [x1_new; x2_new];
end

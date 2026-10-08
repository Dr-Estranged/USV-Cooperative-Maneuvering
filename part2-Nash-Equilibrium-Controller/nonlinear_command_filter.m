function [alpha_f, alpha_dot_f, next_state] = nonlinear_command_filter(alpha, current_state, params, dt)
    % NONLINEAR_COMMAND_FILTER (Levant's Differentiator)
    % 对应原文公式:
    %   dot_alpha_f = -zeta1 * sig(alpha_f - alpha)^0.5 + alpha_d_f
    %   dot_alpha_d_f = -zeta2 * sig(alpha_f - alpha)^0
    
    % 状态解包
    x1 = current_state(1:3); % alpha_f (估计值)
    x2 = current_state(4:6); % alpha_d_f (辅助变量)
    
    % 参数
    z1 = params.zeta1; 
    z2 = params.zeta2; 
    
    % 误差 e = alpha_f - alpha
    err = x1 - alpha;
    
    % [修改] 引入边界层 phi，用平滑函数代替 sign() 以消除高频震荡
    phi = 0.5; % 边界层厚度，调大一点以进一步消除震荡
    soft_sign = tanh(err ./ phi); 
    
    % 计算 sig 函数
    % sig(x)^0.5 = |x|^0.5 * sgn(x) -> 改用软开关
    sig_05 = abs(err).^0.5 .* soft_sign;
    % sig(x)^0 = sgn(x) -> 改用软开关
    sig_0  = soft_sign;
    
    % 计算导数 (Dynamics)
    dx1 = -z1 .* sig_05 + x2;
    dx2 = -z2 .* sig_0;
    
    % 欧拉积分
    x1_new = x1 + dx1 * dt;
    x2_new = x2 + dx2 * dt;
    
    % 输出
    alpha_f = x1_new;
    alpha_dot_f = dx1; % dx1 是 x1 的变化率，即 alpha_dot_f
    
    next_state = [x1_new; x2_new];
end


function dxdt = isv_plant(t, X, tau, model_params)
    % ISV_PLANT 无人艇物理模型
    % 对应论文公式 (1), (2)
    %
    % 运动学方程: eta_dot = R(psi) * nu
    % 动力学方程: M * nu_dot = tau + f(nu) + w(t)
    %
    % 输入:
    %   t: 当前时间
    %   X: [6x1] 状态向量 [x; y; psi; u; v; r]
    %   tau: [3x1] 控制输入 [tau_u; tau_v; tau_r]
    %   model_params: 模型参数结构体
    %
    % 输出:
    %   dxdt: [6x1] 状态导数
    
    %% 状态提取
    psi = X(3);       % 艏向角
    nu  = X(4:6);     % 体坐标系速度 [u; v; r]
    
    % 获取质量矩阵逆
    M_inv = model_params.M_inv;

    %% 1. 运动学方程 (公式 1)
    % eta_dot = R(psi) * nu
    
    R = [cos(psi), -sin(psi), 0;
         sin(psi),  cos(psi), 0;
         0,         0,        1];
     
    eta_dot = R * nu;

    %% 2. 动力学方程 (公式 2)
    % M * nu_dot = tau + f(nu) + w(t)
    
    % --- f(nu): 未知非线性项 (水动力阻力等) ---
    % 线性阻力
    D_lin = diag([50, 50, 10]); 
    % 非线性阻力 (二次阻力)
    D_nonlin = diag([10*abs(nu(1)), 10*abs(nu(2)), 5*abs(nu(3))]);
    
    % 总阻力 (与运动方向相反)
    f_nu = -(D_lin + D_nonlin) * nu; 
    
    % --- w(t): 环境干扰 ---
    w_t = [10 * sin(0.1*t);   % 纵向干扰
           10 * cos(0.1*t);   % 横向干扰
           1 * sin(0.05*t)];  % 转艏干扰
       
    % --- 加速度计算 ---
    nu_dot = M_inv * (tau + f_nu + w_t);

    %% 3. 输出
    dxdt = [eta_dot; nu_dot];
end

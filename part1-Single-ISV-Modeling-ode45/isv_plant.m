function dxdt = isv_plant(t, X, tau, model_params)
    % ISV_PLANT 模拟无人艇物理模型 (对应截图中的公式)
    %
    % 输入 Inputs:
    %   t:            当前时间 (scalar)
    %   X:            当前状态向量 [6x1], 定义为 [x; y; psi; u; v; r]
    %                 其中 eta = [x; y; psi], nu = [u; v; r]
    %   tau:          控制输入 [3x1], 定义为 [tau_u; tau_v; tau_r] (图中标绿部分)
    %   model_params: 包含模型参数的结构体 (M矩阵, 阻力系数等)
    %
    % 输出 Outputs:
    %   dxdt:         状态导数 [6x1], 定义为 [x_dot; y_dot; psi_dot; u_dot; v_dot; r_dot]

    %% 0. 状态解包 (Unpacking)
    psi = X(3);       % 艏向角 (用于计算旋转矩阵)
    nu  = X(4:6);     % 船体坐标系速度 [u; v; r]
    
    % 获取模型参数
    M_inv = model_params.M_inv; % 质量矩阵的逆 (预先计算好以提高效率)

    %% 1. 计算运动学导数 (对应截图上半部分方程)
    % 公式: eta_dot = R(psi) * nu
    
    % 构建旋转矩阵 R(psi) 
    R = [cos(psi), -sin(psi), 0;
         sin(psi),  cos(psi), 0;
         0,         0,        1];
     
    % 计算位置变化率
    eta_dot = R * nu;

    %% 2. 计算动力学导数 (对应截图下半部分方程)
    % 原公式: M * nu_dot = tau + f(nu) + w(t)
    % 移项后: nu_dot = M_inv * (tau + f(nu) + w(t))
    
    % --- A. 计算 f(nu): 未知非线性项 (在这里作为上帝视角的真值) ---
    % 论文中 f_i(nu) 包含水动力阻力、科里奥利力等 [cite: 112]
    % 这里我们需要模拟一个具体的 f(nu)，通常使用 Fossen 模型:
    
    % 线性阻力矩阵 (示例)
    D_lin = diag([50, 50, 10]); 
    % 非线性阻力矩阵 (简单的二次阻力)
    D_nonlin = diag([10*abs(nu(1)), 10*abs(nu(2)), 5*abs(nu(3))]);
    
    % 总阻力 f(nu) (注意符号，通常阻力与运动方向相反)
    f_nu = -(D_lin + D_nonlin) * nu; 
    
    % --- B. 计算 w(t): 外部环境干扰 [cite: 112] ---
    % 假设有一个随时间变化的正弦波干扰
    w_t = [10 * sin(0.1*t); %此时干扰力慢慢增加，船的速度v迅速调整，使得阻力 $f(\nu)$ 几乎实时地抵消了干扰力。
        %如果把0.1改成更高的频率，干扰变化得比船的反应速度快，加速度就会显现出来，变成剧烈震荡。
           10 * cos(0.1*t); 
           1 * sin(0.05*t)];
       
    % --- C. 计算加速度 nu_dot ---
    % 核心公式实现
    nu_dot = M_inv * (tau + f_nu + w_t);

    %% 3. 打包输出
    dxdt = [eta_dot; nu_dot];
end
function dxdt = isv_plant(t, X, tau, model_params)
    % ISV_PLANT 模拟无人艇物理模型 (包含4种独立扰动测试)
    
    % 0. 状态解包
    psi = X(3);       
    nu  = X(4:6);     
    u = nu(1); v = nu(2); r = nu(3);
    
    M_inv = model_params.M_inv; 

    % 1. 计算运动学导数
    R = [cos(psi), -sin(psi), 0;
         sin(psi),  cos(psi), 0;
         0,         0,        1];
    eta_dot = R * nu;

    % 2. 注入外部扰动 w(t) 和内部未建模动态 f(nu)
    w_t = zeros(3, 1);
    f_nu = zeros(3, 1);
    
    % 根据传入的 dist_type 切换扰动工况
    switch model_params.dist_type
        case 1 % (1) 常数型扰动 (Constant)
            w_t = [15; 10; 2];
            
        case 2 % (2) 周期型扰动 (Periodic)
            w_t = [15 * sin(0.2 * t); 
                   15 * cos(0.2 * t); 
                    3 * sin(0.1 * t)];
                    
        case 3 % (3) 随机扰动 (Random Gaussian Noise)
            % randn 生成标准正态分布 N(0,1)
            w_t = [5 * randn(); 
                   5 * randn(); 
                   1 * randn()];
                   
        case 4 % (4) 状态依赖型扰动 (State-dependent)
            % 包含线性阻尼和非线性二次阻尼
            f_nu = [-20*u - 10*abs(u)*u; 
                    -20*v - 10*abs(v)*v; 
                     -5*r -  2*abs(r)*r];
    end

    % 3. 计算动力学导数
    nu_dot = M_inv * (tau + f_nu + w_t);

    % 4. 打包输出
    dxdt = [eta_dot; nu_dot];
end
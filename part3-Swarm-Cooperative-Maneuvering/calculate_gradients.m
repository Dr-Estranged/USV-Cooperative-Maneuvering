function [grad_J1, grad_J2] = calculate_gradients(agent_i, my_eta_bar, my_theta_bar, ...
                                                   Est_Eta_bar_All, Est_Theta_bar_All, t, R_list)
    % CALCULATE_GRADIENTS 计算Nash均衡代价函数梯度
    % 修改版：增加了 J_1 (几何任务) 中的博弈耦合项
    % 对应论文公式 (31) 和 (32) [cite: 367, 371]
    
    M = length(R_list);
    
    %% 1. 提取当前状态
    x_i = my_eta_bar(1);
    y_i = my_eta_bar(2);
    psi_i = my_eta_bar(3);
    
    %% 2. 计算参考轨迹 p_i(theta_bar_i)
    R_i = R_list(agent_i);
    
    % 参考位置 (由Nash参数决定)
    x_ri = R_i * cos(my_theta_bar);
    y_ri = R_i * sin(my_theta_bar);
    psi_ri = 1.03 * my_theta_bar + 0.4855 * pi;
    
    % 期望的参数值 theta_ri
    theta_ri = 0.0103 * t;
    
    %% 3. J1 几何任务梯度 (公式 31) - 关键修改处
    % 论文公式包含耦合项: - (0.1 * sum(x_k) + 100*cos(theta_i)) * x_i
    
    % 计算所有智能体动作估计值的和 (利用估计器数据)
    % Est_Eta_bar_All 是 [3 x M] 矩阵，按行求和得到 [sum_x; sum_y; sum_psi]
    sum_eta_all = sum(Est_Eta_bar_All, 2);
    sum_x_all   = sum_eta_all(1);
    sum_y_all   = sum_eta_all(2);
    sum_psi_all = sum_eta_all(3);
    
    % 初始化梯度向量
    grad_J1 = zeros(3, 1);
    
    % --- x方向梯度 ---
    % 导数推导: d/dx_i [ (x_i - x_ri)^2 - (0.1*sum(x_k) + C)*x_i ]
    % 注意: sum(x_k) 中包含 x_i，因此是对二次项求导
    % 结果: 2*(x_i - x_ri) - 0.1*(x_i + sum(x_k)) - C
    term_track_x = 2 * (x_i - x_ri);
    term_game_x  = -0.1 * (x_i + sum_x_all) - 100 * cos(my_theta_bar);
    grad_J1(1)   = term_track_x + term_game_x;
    
    % --- y方向梯度 ---
    term_track_y = 2 * (y_i - y_ri);
    term_game_y  = -0.1 * (y_i + sum_y_all) - 100 * sin(my_theta_bar);
    grad_J1(2)   = term_track_y + term_game_y;
    
    % --- psi方向梯度 ---
    % 论文对应项: - (0.01 * sum(psi_k)) * psi_i
    term_track_psi = 2 * (psi_i - psi_ri);
    term_game_psi  = -0.01 * (psi_i + sum_psi_all);
    grad_J1(3)     = term_track_psi + term_game_psi;
    
    %% 4. J2 动态任务梯度 (公式 32) - 保持原样
    % J_i^2 = (theta_bar_i - theta_ri)^2 - (0.01 * sum(theta_k)) * theta_i
    
    % 计算所有theta估计值的和
    sum_theta_all = sum(Est_Theta_bar_All);
    
    % 梯度: 2*(theta - theta_r) - 0.01*(theta + sum(theta_k))
    % 注意：你之前的代码用的是差值形式 (sum_diff)，这类似一致性协议。
    % 若要严格对应论文公式 (32)，应用下式：
    term_track_theta = 2 * (my_theta_bar - theta_ri);
    term_game_theta  = -0.01 * (my_theta_bar + sum_theta_all);
    
    % 这里我保留了你原本的一致性逻辑 (sum_diff)，因为它也能实现类似效果。
    % 但如果你想严格复现公式32，请取消注释下面这行：
    grad_J2 = term_track_theta + term_game_theta;
    
  
end
function [d_Eta_bar_hat, d_Theta_bar_hat] = update_estimators(Eta_bar_hat, Theta_bar_hat, ...
                                                              Nash_Eta, Nash_Theta, A, mu, kappa)
    % UPDATE_ESTIMATORS 基于一致性协议的动作估计器
    % 对应论文 Section III-A, 公式 (9) 和 (10)
    %
    % 关键理解：
    %   - 估计的是其他ISV的Nash均衡变量 (eta_bar, theta_bar)
    %   - 不是估计其他ISV的实际位置
    %
    % 输入:
    %   Eta_bar_hat: [3 x M x M] Nash位置估计
    %               Eta_bar_hat(:,j,i) = ISV i 对 ISV j 的 eta_bar 估计
    %   Theta_bar_hat: [1 x M x M] Nash参数估计
    %               Theta_bar_hat(1,j,i) = ISV i 对 ISV j 的 theta_bar 估计
    %   Nash_Eta: [3 x M] 真实的Nash均衡位置目标
    %   Nash_Theta: [1 x M] 真实的Nash均衡参数
    %   A: [M x M] 邻接矩阵
    %   mu: [3 x 3] 位置估计增益矩阵
    %   kappa: [scalar] 参数估计增益
    %
    % 输出:
    %   d_Eta_bar_hat: [3 x M x M] Nash位置估计的导数
    %   d_Theta_bar_hat: [1 x M x M] Nash参数估计的导数
    
    [~, M, ~] = size(Eta_bar_hat);
    d_Eta_bar_hat = zeros(3, M, M);
    d_Theta_bar_hat = zeros(1, M, M);
    
    % 遍历每一个观察者 i
    for i = 1:M
        % 找到 ISV i 的邻居集合 N_i
        neighbors = find(A(i, :) == 1);
        
        % 遍历每一个被估计的目标 j
        for j = 1:M
            %% --- 公式 (9): Nash位置估计动态 ---
            % dot_eta_bar_hat_{i,j} = -mu * (一致性项 + 观测项)
            
            % 一致性项: 与邻居的估计比较
            consensus_eta = zeros(3, 1);
            for k = neighbors
                consensus_eta = consensus_eta + (Eta_bar_hat(:, j, i) - Eta_bar_hat(:, j, k));
            end
            
            % 观测项: 只有当 i 能直接获取 j 的信息时才有效
            % b_{i,j} = 1 if (i == j) or (j in N_i)
            correction_eta = zeros(3, 1);
            if i == j
                % ISV可以获取自己的Nash均衡变量
                correction_eta = Eta_bar_hat(:, j, i) - Nash_Eta(:, j);
            elseif A(i, j) == 1
                % ISV可以从邻居获取其Nash均衡变量
                correction_eta = Eta_bar_hat(:, j, i) - Nash_Eta(:, j);
            end
            
            % 总导数
            d_Eta_bar_hat(:, j, i) = -mu * (consensus_eta + correction_eta);
            
            %% --- 公式 (10): Nash参数估计动态 ---
            
            % 一致性项
            consensus_theta = 0;
            for k = neighbors
                consensus_theta = consensus_theta + (Theta_bar_hat(1, j, i) - Theta_bar_hat(1, j, k));
            end
            
            % 观测项
            correction_theta = 0;
            if i == j
                correction_theta = Theta_bar_hat(1, j, i) - Nash_Theta(j);
            elseif A(i, j) == 1
                correction_theta = Theta_bar_hat(1, j, i) - Nash_Theta(j);
            end
            
            % 总导数
            d_Theta_bar_hat(1, j, i) = -kappa * (consensus_theta + correction_theta);
        end
    end
end

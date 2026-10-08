function dual_gradient_bar3()
    % 1. 准备数据
    [X, Y] = meshgrid(1:15, 1:12);
    % 制造一些有起伏的数据
    Z = 2000 + X.*100 + Y.*50 + rand(12,15)*500; 
    % 让右半部分稍微高一点，区分度更明显
    Z(:, 10:end) = Z(:, 10:end) + 1000;

    % 2. 绘制基础 3D 柱状图
    figure('Color', 'w', 'Position', [100, 100, 800, 600]);
    h = bar3(Z, 0.6); % 0.6 是柱子宽度
    
    % ==========================================================
    % 核心部分：制作双色渐变 Colormap
    % ==========================================================
    n_levels = 128; % 渐变级数，越高越细腻
    
    % --- 定义青色渐变 (由浅白 -> 深青) ---
    % RGB: 浅[0.9, 1, 1] -> 深[0.1, 0.6, 0.7]
    cyan_map = [linspace(0.9, 0.1, n_levels)', ...
                linspace(1.0, 0.6, n_levels)', ...
                linspace(1.0, 0.7, n_levels)'];
            
    % --- 定义粉色渐变 (由浅白 -> 深红) ---
    % RGB: 浅[1, 0.9, 0.9] -> 深[0.8, 0.3, 0.3]
    pink_map = [linspace(1.0, 0.8, n_levels)', ...
                linspace(0.9, 0.3, n_levels)', ...
                linspace(0.9, 0.3, n_levels)'];

    % --- 拼接色卡 ---
    full_map = [cyan_map; pink_map];
    colormap(full_map); % 应用这个 256 色的长色卡
    
    % ==========================================================
    % 核心部分：修改每个柱子的颜色映射 (CData)
    % ==========================================================
    
    % 找到数据的最大高度，用于归一化颜色，保证最高的柱子最深
    max_z = max(Z(:)); 
    min_z = 0; % 假设柱子从0开始
    
    split_index = 9; % 在第9列分界
    
    for i = 1:length(h)
        % 获取当前列的 Z 坐标数据 (用于计算颜色深浅)
        zData = h(i).ZData; 
        
        % 归一化：将高度映射到 0 到 1 之间
        % 这样底座就是 0 (浅色)，顶端就是 1 (深色)
        cData = (zData - min_z) / (max_z - min_z);
        
        % 关键：根据列号 i，决定去色卡的哪一部分取色
        if i <= split_index
            % --- 左侧：使用色卡的前半段 (1 ~ 128) ---
            % 映射公式：Index = 1 + cData * (n_levels - 1)
            h(i).CData = 1 + cData * (n_levels - 1);
        else
            % --- 右侧：使用色卡的后半段 (129 ~ 256) ---
            % 映射公式：Index = (n_levels + 1) + cData * (n_levels - 1)
            h(i).CData = (n_levels + 1) + cData * (n_levels - 1);
        end
        
        % --- 开启渐变渲染 ---
        h(i).CDataMapping = 'direct'; % 告诉 MATLAB 直接使用我们计算的索引
        h(i).FaceColor = 'interp';    % 开启插值渲染 (实现平滑渐变)
        h(i).EdgeColor = [0.6 0.6 0.6]; % 边框设为淡灰色，避免干扰渐变
        h(i).FaceAlpha = 0.9; % 稍微透明一点点，质感更好
    end
    
    % ==========================================================
    % 装饰图表
    % ==========================================================
    view(-40, 25);
    zlabel('Height / m');
    ylabel('Time / s');
    xlabel('Velocity / m/s');
    grid on;
    
    % 修正坐标轴刻度文本(示例)
    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.YDir = 'reverse';
end
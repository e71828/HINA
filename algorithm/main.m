%% 计算子载波间隔
% 参数定义
fC = 2565e6; % 载频，单位Hz
B = 100e6; % 带宽，单位Hz
scs = 30e3; % 子载波间隔，单位Hz
comb_spacing = 4; % comb间隔，每4个子载波放置1个SRS
num_srs_subcarriers = 816; % SRS的有效子载波数
TC = 1/(480 * 1000 * 4096); 

N = B / scs;                % 子载波数量
delta_f = B / N;            % 子载波带宽
srs_spacing = comb_spacing * scs;  % SRS信号的频率间隔
%% 读取输入文件
pilot = load("../pilot and example/pilot.mat");
Xf = pilot.pilot;
%% 读取文件夹
folder = '../data';  % 文件夹路径
filePattern = fullfile(folder, '*.mat');  % 指定文件类型，这里是MAT文件
matFiles = dir(filePattern);  % 获取所有符合要求的文件信息
tau_ans = zeros(1,800);
%% 处理前400个文件
for i = 1:400
    filename = fullfile(matFiles(i).folder, matFiles(i).name);  % 获取文件名及路径
    data = load(filename);  % 加载MAT文件中的数据
    variable_names = who('-file', filename);  % 获取MAT文件中的变量名
    variable_name = variable_names{1};  % 假设MAT文件中只有一个变量
    Yf = data.(variable_name);  % 获取MAT文件中的变量值    % 对数据进行处理
    x = Yf./Xf;
    X = abs(ifft(x));
    Nx = length(x); % 信号长度
    X_sel = X(1:256);
    t_sel = (0:255)*(1/srs_spacing/Nx)/TC;

    % 使用findpeaks函数找到所有峰值
    [peaks, locs] = findpeaks(X_sel);
    
    % 获取最大峰值
    max_peak = max(peaks);
    
    % 筛选出幅度大于等于最大峰值一半的峰值，并且不在第一个值的位置
    valid_peaks = peaks(peaks >= max_peak / 3 & locs ~= 1);
    valid_locs = locs(peaks >= max_peak / 3 & locs ~= 1);
    
    % 如果存在有效峰值，获取第一个有效峰值及其位置
    if ~isempty(valid_peaks)
        first_peak = valid_peaks(1);
        first_loc = valid_locs(1);
        
        % 获取峰值左侧、峰值、峰值右侧的数据
        peak_data = X_sel(first_loc-1:first_loc+1);
        time_data = t_sel(first_loc-1:first_loc+1);
    
        % 使用二次多项式拟合峰值附近的数据
        p = polyfit(time_data, peak_data, 2);
    
        % 计算二次多项式的顶点（即真实峰值位置）
        estimated_t = -p(2) / (2 * p(1));
    
        % 输出结果
        % fprintf('Estimated true time of first valid peak: %f\n', estimated_t);
    else
        fprintf('No valid peaks found.\n');
    end
    tau_ans(i) = estimated_t;
end 
isInRange = (tau_ans(1:400) >= 0) & (tau_ans(1:400) <= 190); % 
allInRange = all(isInRange, 'all');
disp(find(~isInRange));
assert(allInRange, 'Not all elements are in the range.');
%% 处理后400个文件
tau_est = zeros(4,1);
for i = 401:800
    filename = fullfile(matFiles(i).folder, matFiles(i).name);  % 获取文件名及路径
    data = load(filename);  % 加载MAT文件中的数据
    variable_names = who('-file', filename);  % 获取MAT文件中的变量名
    variable_name = variable_names{1};  % 假设MAT文件中只有一个变量
    Yf = data.(variable_name);  % 获取MAT文件中的变量值    % 对数据进行处理
    x = Yf./Xf;
    X = abs(ifft(x,[],2));
    Nx = length(x); % 信号长度
    X_sel_4 = X(:,1:256);
    t_sel = (0:255)*(1/srs_spacing/Nx)/TC;
    for j = 1:4
        X_sel = X_sel_4(j,:);
    
        % 使用findpeaks函数找到所有峰值
        [peaks, locs] = findpeaks(X_sel);
        
        % 获取最大峰值
        max_peak = max(peaks);
        
        % 筛选出幅度大于等于最大峰值一半的峰值，并且不在第一个值的位置
        valid_peaks = peaks(peaks >= max_peak / 3 & locs ~= 1);
        valid_locs = locs(peaks >= max_peak / 3 & locs ~= 1);
        
        % 如果存在有效峰值，获取第一个有效峰值及其位置
        if ~isempty(valid_peaks)
            first_peak = valid_peaks(1);
            first_loc = valid_locs(1);
            
            % 获取峰值左侧、峰值、峰值右侧的数据
            peak_data = X_sel(first_loc-1:first_loc+1);
            time_data = t_sel(first_loc-1:first_loc+1);
        
            % 使用二次多项式拟合峰值附近的数据
            p = polyfit(time_data, peak_data, 2);
        
            % 计算二次多项式的顶点（即真实峰值位置）
            estimated_t = -p(2) / (2 * p(1));
        
            % 输出结果
            % fprintf('Estimated true time of first valid peak: %f\n', estimated_t);
        else
            fprintf('No valid peaks found.\n');
        end
        tau_est(j) = estimated_t;
    end
    tau_ans(i) = tau_est_merge(tau_est);
end 
%% 写入答案文件
fileID = fopen('../data/answer.txt', 'w');
for value = tau_ans
    fprintf(fileID, "%.2f,\r\n", value);
end
fclose(fileID);
%% 合并文件

%% 压缩文件
cd ..
zip('xxx-西北工业大学.zip', {'algorithm', 'data/answer.txt'});
unzip('xxx-西北工业大学.zip','xxx-西北工业大学');
delete('xxx-西北工业大学.zip');
%% 四个天线估计融合
function result = tau_est_merge(tau_est)
    % 假设tau_est是一个4x1的矩阵，包含四个估计值
    
    % 计算标准差
    standard_deviation = std(tau_est);
    
    % 初始化结果变量
    result = NaN;
    
    % 检查标准差是否小于6
    if standard_deviation < 6
        % 如果标准差小于6，则取均值
        result = mean(tau_est);
        
    % 检查标准差是否大于6且小于20
    elseif standard_deviation > 6 && standard_deviation < 20
        % 如果标准差大于6且小于20，则从最大到最小划分3段，选择两个数的数据段，取其中两数的中值
        
        % 对tau_est进行排序
        sorted_tau_est = sort(tau_est, 'descend');
        
        % 划分3段
        num_segments = 3;
        segment_size = ceil(length(sorted_tau_est) / num_segments);
        
        % 选择两个数的数据段
        two_number_segment = sorted_tau_est(segment_size + 1 : 2 * segment_size);
        
        % 计算两个数的中值
        result = median(two_number_segment);
        
    % 检查标准差是否大于20
    elseif standard_deviation > 20
        % 如果标准差大于20，则取最小两个值的中值
        
        % 对tau_est进行排序
        sorted_tau_est = sort(tau_est);
        
        % 取最小两个值
        min_two_values = sorted_tau_est(1:2);
        
        % 计算最小两个值的中值
        result = median(min_two_values);
    end
end

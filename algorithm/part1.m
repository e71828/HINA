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
for i = 1:20
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
    % 绘制谱
    figure;
    plot(t_sel,abs(X_sel));
    grid on;
    xlim([0 256])

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
        fprintf('Estimated true time of first valid peak: %f\n', estimated_t);
    else
        fprintf('No valid peaks found.\n');
    end
    tau_ans(i) = estimated_t;
end 
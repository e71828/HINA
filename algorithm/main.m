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
tau_ans2 = zeros(4,800);

M = 250;       % 协方差矩阵的阶数
N_fft = 32768; % FFT点数（用于计算谱估计）
%% 处理前400个文件
for i = 1:400
    filename = fullfile(matFiles(i).folder, matFiles(i).name);  % 获取文件名及路径
    data = load(filename);  % 加载MAT文件中的数据
    variable_names = who('-file', filename);  % 获取MAT文件中的变量名
    variable_name = variable_names{1};  % 假设MAT文件中只有一个变量
    Yf = data.(variable_name);  % 获取MAT文件中的变量值    % 对数据进行处理
    Hf = Yf./Xf;
    Nsig = mdltest_mcov(Hf');

    % 调用MUSIC算法进行谱估计（不绘制谱估计结果）
    [f_est, P_music] = music_algorithm(Hf, M, Nsig, N_fft);
    
    % 延迟为正，频率为负，反转谱序列
    P_music = P_music(end:-1:1);
    
    % 寻找峰值
    [~, peak_indices] = findpeaks(P_music, 'SortStr', 'descend', 'NPeaks', Nsig);
    f_est_peaks = f_est(peak_indices);
    tau_ans(i) =  f_est_peaks(1)/TC/srs_spacing;
end 
%% 处理后400个文件
tau_est = zeros(4,1);
for i = 401:800
    filename = fullfile(matFiles(i).folder, matFiles(i).name);  % 获取文件名及路径
    data = load(filename);  % 加载MAT文件中的数据
    variable_names = who('-file', filename);  % 获取MAT文件中的变量名
    variable_name = variable_names{1};  % 假设MAT文件中只有一个变量
    Yf = data.(variable_name);  % 获取MAT文件中的变量值    % 对数据进行处理
    Hf = Yf./Xf;
    for j = 1:size(Hf,1)
        Nsig = mdltest_mcov(Hf(j,:)');

        % 调用MUSIC算法进行谱估计（不绘制谱估计结果）
        [f_est, P_music] = music_algorithm(Hf(j,:), M, Nsig, N_fft);
        
        % 延迟为正，频率为负，反转谱序列
        P_music = P_music(end:-1:1);
        
        % 寻找峰值
        [~, peak_indices] = findpeaks(P_music, 'SortStr', 'descend', 'NPeaks', Nsig);
        f_est_peaks = f_est(peak_indices);
        tau_est(j) = f_est_peaks(1)/TC/srs_spacing;
    end
    tau_ans2(:,i) = tau_est';

    % 计算 tau_est 的标准差
    tau_std = std(tau_est);
    
    % 根据 tau_est 的标准差计算 tau_ans 的值
    if tau_std < 10
        % 如果标准差小于 10，则计算 tau_est 的均值
        tau_ans(i) = mean(tau_est);
    else
        % 否则，取 tau_est 的最小值
        tau_ans(i) = min(tau_est);
        fprintf('the std is too large of No. %d\n', i);
    end
end
%% 写入答案文件
if isfile('tau_ans1.mat')
    load tau_ans1.mat
else
    save tau_ans2 tau_ans2
    save tau_ans1 tau_ans1
end
fileID = fopen('../data/answer.txt', 'w');
for value = tau_ans
    fprintf(fileID, "%.2f,\r\n", value);
end
fclose(fileID);
%% 压缩文件
cd ..
zip('xxx-西北工业大学.zip', {'algorithm', 'data/answer.txt'});
unzip('xxx-西北工业大学.zip','唐金峰-西北工业大学');
delete('xxx-西北工业大学.zip');
zip('唐金峰-西北工业大学.zip', '唐金峰-西北工业大学');
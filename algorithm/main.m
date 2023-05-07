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
tau_ans2 = zeros(4,400);

N_fft = 32768; % FFT点数（用于计算谱估计）
%% 处理前400个文件
parfor i = 1:400
    filename = fullfile(matFiles(i).folder, matFiles(i).name);  % 获取文件名及路径
    data = load(filename);  % 加载MAT文件中的数据
    % variable_names = who('-file', filename);  % 获取MAT文件中的变量名
    % variable_name = variable_names{1};  % 假设MAT文件中只有一个变量
    variable_name = 'ant1_data';
    Yf = data.(variable_name);  % 获取MAT文件中的变量值    % 对数据进行处理
    Hf = Yf./Xf;
    [Nsig, R] = mdltest_mcov(Hf', 'fb');
    % 调用MUSIC算法进行谱估计（不绘制谱估计结果）
    [f_est, P_music] = music_algorithm(R, Nsig, N_fft, false);
    % 寻找峰值
    [peak_values, peak_indices] = findpeaks(P_music, 'SortStr', 'descend', 'NPeaks', Nsig);
    peak_indices(peak_values<max(peak_values)-20) = [];
    f_est_peaks = f_est(peak_indices);
    tau_ans(i) =  min(f_est_peaks)/TC/srs_spacing;
end 
%% 处理后400个文件
parfor i = 401:800
    filename = fullfile(matFiles(i).folder, matFiles(i).name);  % 获取文件名及路径
    data = load(filename);  % 加载MAT文件中的数据
    % variable_names = who('-file', filename);  % 获取MAT文件中的变量名
    % variable_name = variable_names{1};  % 假设MAT文件中只有一个变量
    variable_name = 'ant4_data';
    Yf = data.(variable_name);  % 获取MAT文件中的变量值    % 对数据进行处理
    Hf = Yf./Xf;
    P_music_set = zeros(4,1024);
    Nsig = zeros(4,1);
    f_est = linspace(0, 1, N_fft);
    f_est = f_est(1:1024);
    for j = 1:size(Hf,1)
        [Nsig(j), R] = mdltest_mcov(Hf(j,:)', 'fb');
        % 调用MUSIC算法进行谱估计（不绘制谱估计结果）
        [~, P_music] = music_algorithm(R, Nsig(j), N_fft, false);
        P_music_set(j,:) = P_music;
    end
    tau_est = zeros(4,1);
    for j = 1:4
        [peak_values, peak_indices] = findpeaks(P_music_set(j,:), 'SortStr', 'descend', 'NPeaks', max(Nsig));
        peak_indices(peak_values<max(peak_values)-20) = [];
        f_est_peaks = f_est(peak_indices);
        tau_est(j) = min(f_est_peaks)/TC/srs_spacing;
    end
    P_music = mean(P_music_set,1);
    [peak_values, peak_indices] = findpeaks(P_music, 'SortStr', 'descend', 'NPeaks', max(Nsig));
    peak_indices(peak_values<max(peak_values)-20) = [];
    f_est_peaks = f_est(peak_indices);
    tau_ans(i) = min(f_est_peaks)/TC/srs_spacing;
    tau_ans2(:,i-400) = tau_est;
end
%% 保存数据
save tau_ans2 tau_ans2
save tau_ans1 tau_ans
%% 写入答案文件
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
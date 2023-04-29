%% 要求 R2018b之后。
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
close all
folder = '../data';  % 文件夹路径
filePattern = fullfile(folder, '*.mat');  % 指定文件类型，这里是MAT文件
matFiles = dir(filePattern);  % 获取所有符合要求的文件信息
tau_ans = zeros(1,800);
tau_ans2 = zeros(4,400);

M = 250;       % 协方差矩阵的阶数
N_fft = 32768; % FFT点数（用于计算谱估计）
%% 处理后400个文件
for i = (14)+400
    f_est = linspace(0, 1, N_fft);
    filename = fullfile(matFiles(i).folder, matFiles(i).name);  % 获取文件名及路径
    data = load(filename);  % 加载MAT文件中的数据
    % variable_names = who('-file', filename);  % 获取MAT文件中的变量名
    % variable_name = variable_names{1};  % 假设MAT文件中只有一个变量
    variable_name = 'ant4_data';
    Yf = data.(variable_name);  % 获取MAT文件中的变量值    % 对数据进行处理
    Hf = Yf./Xf;
    P_music_set = zeros(4,N_fft);
    Nsig = zeros(4,1);
    for j = 1:size(Hf,1)
        Nsig(j) = mdltest_mcov(Hf(j,:)');
        % 调用MUSIC算法进行谱估计（不绘制谱估计结果）
        [~, P_music] = music_algorithm(Hf(j,:), M, Nsig(j), N_fft);
        % 延迟为正，频率为负，反转谱序列
        P_music = P_music(end:-1:1);
        P_music_set(j,:) = P_music;
    end
    %% 画图
    tau_est = zeros(4,1);
    for j = 1:4
        [~, peak_indices] = findpeaks(P_music_set(j,:), 'SortStr', 'descend', 'NPeaks', min(Nsig));
        f_est_peaks = f_est(peak_indices);
        % 输出估计的频率和真实频率
        disp('Estimated Frequencies:');
        disp(f_est_peaks/TC/srs_spacing);
    
        % 绘制MUSIC谱估计结果
        figure;
        plot(f_est/TC/srs_spacing, P_music_set(j,:), 'LineWidth', 1.2);
        xlabel('Frequency (units of pi)');
        ylabel('Magnitude / dB');
        title('MUSIC Spectrum');
        grid on;
    
        % 标记估计的频率位置
        hold on;
        stem(f_est_peaks/TC/srs_spacing, max(P_music_set(j,:)) * ones(size(f_est_peaks)), 'r', 'filled');
        hold off;
        xlim([0 256])
        tau_est(j) = min(f_est_peaks)/TC/srs_spacing;
    end
    P_music = mean(P_music_set,1);
    [~, peak_indices] = findpeaks(P_music, 'SortStr', 'descend', 'NPeaks', min(Nsig));
    f_est_peaks = f_est(peak_indices);
    tau_ans(i) = min(f_est_peaks)/TC/srs_spacing;
    tau_ans2(:,i-400) = tau_est;
end
%% 判断合理
isInRange = (tau_ans(401:800) >= 0) & (tau_ans(401:800) <= 190); % 
allInRange = all(isInRange, 'all');
disp(find(~isInRange)+400);
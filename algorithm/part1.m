%% 计算子载波间隔
% 参数定义
fC = 2565e6; % 载频，单位Hz
B = 100e6; % 带宽，单位Hz
scs = 30e3; % 子载波间隔，单位Hz
comb_spacing = 4; % comb间隔，每4个子载波放置1个SRS
num_srs_subcarriers = 816; % SRS的有效子载波数
TC = 1/(480 * 1000 * 4096);

N = B / scs;                % 子载波数量
delta_f = B / N;            % 子载波带宽w
srs_spacing = comb_spacing * scs;  % SRS信号的频率间隔
%% 读取输入文件
pilot = load("../pilot and example/pilot.mat");
Xf = pilot.pilot;
%% 读取文件夹
folder = '../data';  % 文件夹路径
filePattern = fullfile(folder, '*.mat');  % 指定文件类型，这里是MAT文件
matFiles = dir(filePattern);  % 获取所有符合要求的文件信息
tau_ans = zeros(1,800);

N_fft = 32768; % FFT点数（用于计算谱估计）
%% 处理前400个文件
for i = 100:110
    filename = fullfile(matFiles(i).folder, matFiles(i).name);  % 获取文件名及路径
    data = load(filename);  % 加载MAT文件中的数据
    % variable_names = who('-file', filename);  % 获取MAT文件中的变量名
    % variable_name = variable_names{1};  % 假设MAT文件中只有一个变量
    variable_name = 'ant1_data';
    Yf = data.(variable_name);  % 获取MAT文件中的变量值    % 对数据进行处理
    Hf = Yf./Xf;
    [Nsig, R] = mdltest_mcov(Hf','fb');

    % 调用MUSIC算法进行谱估计（不绘制谱估计结果）
    [f_est, P_music] = music_algorithm(R, Nsig, N_fft, false);

    % 寻找峰值
    [peak_values, peak_indices] = findpeaks(P_music, 'SortStr', 'descend', 'NPeaks', Nsig);
    peak_indices(peak_values<max(peak_values)-10) = [];
    f_est_peaks = f_est(peak_indices);
    % 输出估计的频率和真实频率
    disp('Estimated Frequencies:');
    disp(f_est_peaks/TC/srs_spacing);

    % 绘制MUSIC谱估计结果
    figure;
    plot(f_est/TC/srs_spacing, P_music, 'LineWidth', 1.2);
    xlabel('Frequency (units of TC)');
    ylabel('Magnitude / dB');
    title('MUSIC Spectrum');
    grid on;

    % 标记估计的频率位置
    hold on;
    stem(f_est_peaks/TC/srs_spacing, max(P_music) * ones(size(f_est_peaks)), 'r', 'filled');
    hold off;
    xlim([0 256])

    tau_ans(i) =  min(f_est_peaks)/TC/srs_spacing;
end 
%% 判断合理
isInRange = (tau_ans(1:400) >= 0) & (tau_ans(1:400) <= 190); % 
allInRange = all(isInRange, 'all');
disp(find(~isInRange));
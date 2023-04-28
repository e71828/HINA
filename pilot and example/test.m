%% 利用参数
% 参数定义
fC = 2565e6; % 载频，单位Hz
B = 100e6; % 带宽，单位Hz
scs = 30e3; % 子载波间隔，单位Hz
comb_spacing = 4; % comb间隔，每4个子载波放置1个SRS
num_srs_subcarriers = 816; % SRS的有效子载波数

N = B / scs;                % 子载波数量
delta_f = B / N;            % 子载波带宽
srs_spacing = comb_spacing * scs;  % SRS信号的频率间隔
srs_bw = num_srs_subcarriers * delta_f;  % SRS信号的带宽

%% 1. 群延迟方法
pilot = load("pilot.mat");
example_64Tc = load("example_64Tc.mat");
Xf = pilot.pilot;
Yf = example_64Tc.example_64Tc;

TC = 1/(480 * 1000 * 4096);
tau = -diff(unwrap(angle(Yf./Xf))) / (2*pi*srs_spacing) / TC;  % 计算群延迟
figure
plot(unwrap(angle(Xf))/pi)
hold on
plot(unwrap(angle(Yf./Xf))/pi)
%% 查看发射相位的包络
data = unwrap(angle(Xf))/pi;

% 使用 findpeaks 函数查找局部极大值（峰值）
[peak_values, peak_indices] = findpeaks(data);

% 使用 findpeaks 函数查找局部极小值（谷值）
[valley_values, valley_indices] = findpeaks(-data);
valley_values = -valley_values;

% 首端和末端处理
if data(2)>data(1)
    valley_values = [data(1) valley_values];
    valley_indices = [1 valley_indices];
else
    peak_values = [data(1) peak_values];
    peak_indices = [1 peak_indices];
end
if data(end)>data(end-1)
    peak_values = [peak_values data(end)];
    peak_indices = [peak_indices length(data)];
else
    valley_values = [valley_values data(end)];
    valley_indices = [valley_indices length(data)];
end

% 使用 griddedInterpolant 创建插值对象计算上包络
F_upper = griddedInterpolant(peak_indices, peak_values, 'spline');
upper_envelope = F_upper(1:length(data));

% 使用 griddedInterpolant 创建插值对象计算下包络
F_lower = griddedInterpolant(valley_indices, valley_values, 'spline');
lower_envelope = F_lower(1:length(data));


% 绘制原始数据、上包络、下包络以及峰值数据
figure
plot(data, 'b');  % 原始数据用蓝色绘制
hold on;
plot(upper_envelope, 'r');  % 上包络用红色绘制
plot(lower_envelope, 'g');  % 下包络用绿色绘制
plot(peak_indices, peak_values, 'o');  % 峰值数据用红色圆圈标记
plot(valley_indices, valley_values, 'o');  % 谷值数据用绿色圆圈标记
hold off;
legend('Original Data', 'Upper Envelope', 'Lower Envelope', 'Peak Values', 'Valley Values');
xlabel('Time');
ylabel('Amplitude');
title('Data with Upper and Lower Envelopes and Peak Values');
%% 2. 逆Fourier方法
% generate x with Yf and Xf.
Hf = Yf./Xf;
h = ifft(Hf);

% 计算轴
Nh = length(h); % 信号长度
Ts = 1/srs_spacing;
t = (0:Nh-1)*(Ts/Nh); % 频率轴

% 绘制谱
figure;
plot(t/TC, abs(h));
grid on;
xlim([0 256])

% 应该是 64 TC，但只能得到60 TC。

%% 创建一个复数信号
% 将时延对应相位当成频率值，频率坐标点当成时间序列。对时间信号做处理。
f_true = [64 256];
Hf = sum(2*exp(-1i * 2 * pi * f_true'*TC*srs_spacing .* (1:num_srs_subcarriers)), 1);
Hf = awgn(Hf,0);
figure
plot(unwrap(angle(Hf))/pi)  % 相位展开，难点
hold on
plot(abs(Hf))
Nsig = mdltest_mcov(Hf');
disp(Nsig)
%% 复数信号案例-测试时延分辨率
% 对复数信号进行逆傅里叶变换
h = ifft(Hf);

% 计算轴
Nh = length(h); % 信号长度
t = (0:Nh-1)*(Ts/Nh); % 频率轴

% 绘制谱
figure;
plot(t/TC, abs(h));
grid on;
xlim([0 256])

% 时延分辨率低为 20TC。

%% MUSIC
% 参数设置
M = 250;       % 协方差矩阵的阶数
L = 2;       % 信号子空间的维度（信号源的数量）
N_fft = 32768; % FFT点数（用于计算谱估计）
f_true = [124 36];
Hf = sum(2*exp(-1i * 2 * pi * f_true'*TC*srs_spacing .* (1:num_srs_subcarriers)), 1);

Nsig = mdltest_mcov(Hf');
% 调用MUSIC算法进行谱估计（不绘制谱估计结果）
[f_est, P_music] = music_algorithm(Hf, M, Nsig, N_fft);

% 延迟为正，频率为负，反转谱序列
P_music = P_music(end:-1:1);

% 寻找峰值
[~, peak_indices] = findpeaks(P_music, 'SortStr', 'descend', 'NPeaks', Nsig);
f_est_peaks = f_est(peak_indices);

% 输出估计的频率和真实频率
disp('Estimated Frequencies:');
disp(f_est_peaks/TC/srs_spacing);
disp('True Frequencies:');
disp(f_true);

% 绘制MUSIC谱估计结果
figure;
plot(f_est, P_music, 'LineWidth', 1.2);
xlabel('Frequency (units of pi)');
ylabel('Magnitude / dB');
title('MUSIC Spectrum');
grid on;

% 标记估计的频率位置
hold on;
stem(f_est_peaks, max(P_music) * ones(size(f_est_peaks)), 'r', 'filled');
hold off;

% 标记真实频率位置
hold on;
stem(f_true*TC*srs_spacing, max(P_music) * ones(size(f_true)), 'g', 'filled');
hold off;

legend('MUSIC Spectrum', 'Estimated Frequencies', 'True Frequencies');
xlim([0 5e-3*pi])
%% MUSIC
function [f_est, P_music] = music_algorithm(x, M, L, N_fft, plot_spectrum)
% MUSIC算法进行谱估计
% 输入：
%   x - 输入信号
%   M - 协方差矩阵的阶数
%   L - 信号子空间的维度（信号源的数量）
%   N_fft - FFT点数（用于计算谱估计）
%   plot_spectrum - 是否绘制MUSIC谱估计结果（可选，默认为false）
% 输出：
%   f_est - 估计的频率（归一化频率，单位：pi）
%   P_music - MUSIC谱估计

% 设置默认参数
if nargin < 5
    plot_spectrum = false;
end

% 计算自相关函数
rxx = xcorr(x, M-1, 'biased');

% 构造Hermitian Toeplitz协方差矩阵
Rxx = toeplitz(rxx(M:end));

% 计算协方差矩阵的特征值和特征向量
[V, D] = eig(Rxx);

% 对特征值进行排序，并获取排序后的索引
[~, I] = sort(diag(D));

% 提取噪声子空间的特征向量
V_noise = V(:, I(1:M-L));

% 初始化MUSIC谱估计
P_music = zeros(1, N_fft);

% 计算MUSIC谱估计
for k = 1:N_fft
    f = (k - 1) / N_fft;  % 归一化频率
    e = exp(-1j * 2 * pi * f * (0:M-1)).';  % 构造复指数向量
    P_music(k) = 1 / real(e' * (V_noise * V_noise') * e);  % MUSIC谱估计(取实部）
end

% 转换为分贝（dB）表示
P_music = pow2db(P_music);

% 估计的频率范围
f_est = linspace(0, 1, N_fft);

% 如果需要绘制MUSIC谱估计结果，则进行绘图
if plot_spectrum
    figure;
    plot(f_est, P_music, 'LineWidth', 1.2);
    xlabel('Frequency (units of pi)');
    ylabel('Magnitude / dB');
    title('MUSIC Spectrum');
    grid on;
end

end
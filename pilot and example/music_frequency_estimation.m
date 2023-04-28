function music_frequency_estimation()
% 参数设置
N = 64;      % 信号长度
M = 30;       % 协方差矩阵的阶数
L = 2;       % 信号子空间的维度（信号源的数量）
N_fft = 1024; % FFT点数（用于计算谱估计）

% 真实频率（归一化）
f_true = [0.1, 0.01];

% 生成多信号
s = sum(2*exp(1i * 2 * pi * f_true' .* (1:N) + 1j*2*pi*rand(1)), 1);

% 添加复数高斯白噪声
x = s + wgn(1, N, 0, 'complex');

% 调用MUSIC算法进行谱估计（不绘制谱估计结果）
[f_est, P_music] = music_algorithm(x, M, L, N_fft);

[~, peak_indices] = findpeaks(P_music, 'SortStr', 'descend', 'NPeaks', L);
f_est_peaks = f_est(peak_indices);

% 输出估计的频率和真实频率
disp('Estimated Frequencies:');
disp(f_est_peaks);
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
stem(f_true, max(P_music) * ones(size(f_true)), 'g', 'filled');
hold off;

legend('MUSIC Spectrum', 'Estimated Frequencies', 'True Frequencies');

end


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

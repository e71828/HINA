%% MUSIC
function [f_est, P_music] = music_algorithm(x, M, L, N_fft, plot_spectrum, fb)
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
    fbFlag = false;
end

if  nargin == 6
    validatestring(fb,{'fb'},'music_algorithm','',6);
    fbFlag = true;
end
% 计算自相关函数
rxx = xcorr(x, M-1, 'biased');

% 构造Hermitian Toeplitz协方差矩阵
Rxx = toeplitz(rxx(M:end));

if fbFlag
  Rxx = spsmooth(Rxx,1,'fb');
end

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
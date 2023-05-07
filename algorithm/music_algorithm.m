%% MUSIC
function [f_est, P_music] = music_algorithm(R, M, N_fft, plot_spectrum)
% MUSIC算法进行谱估计
% 输入：
%   R - 输入信号的自相关
%   M - 信号子空间的维度（信号源的数量）
%   N_fft - FFT点数（用于计算谱估计）
%   plot_spectrum - 是否绘制MUSIC谱估计结果（可选，默认为false）
% 输出：
%   f_est - 估计的频率（归一化频率，单位：pi）
%   P_music - MUSIC谱估计

% 设置默认参数
if nargin < 4
    plot_spectrum = false;
end


Rxx = R;
N = length(R);

% 计算协方差矩阵的特征值和特征向量
[V, D] = eig(Rxx);

% 对特征值进行排序，并获取排序后的索引
[~, I] = sort(diag(D));


% 提取噪声子空间的特征向量
V_noise = V(:, I(1:N-M));

% 初始化MUSIC谱估计
P_music = zeros(1, N_fft);

% 估计的频率范围
f_est = linspace(0, 1, N_fft);

% 计算MUSIC谱估计
for k = N_fft-1023:N_fft
    f = f_est(k);  % 归一化频率
    e = exp(-1j * 2 * pi * f * (0:N-1)).';  % 构造复指数向量
    P_music(k) = 1 / real(e' * (V_noise * V_noise') * e);  % MUSIC谱估计(取实部）
end

% 延迟为正，频率为负，反转谱序列
P_music = P_music(end:-1:1);
P_music = P_music(1:1024);
f_est = f_est(1:1024);

% 转换为分贝（dB）表示
P_music = pow2db(P_music);



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
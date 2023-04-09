%% 1. 群延迟方法
pilot = load("pilot.mat");
example_64Tc = load("example_64Tc.mat");
Xf = pilot.pilot;
Yf = example_64Tc.example_64Tc;

% SRS的频率间隔由下一 section 求得。
TC = 1/(480 * 1000 * 4096); 
tau = -diff(unwrap(angle(Yf./Xf))) / (2*pi*srs_spacing) / TC;  % 计算群延迟

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

% 计算频率中心，估计值
srs_subcarriers = floor((N - num_srs_subcarriers) / 2) + (1:num_srs_subcarriers);
frequencies = fC + (srs_subcarriers - N/2) * delta_f;

%% 2. 逆Fourier方法
% generate x with Yf and Xf.
x = Yf./Xf;
X = ifft(x);

% 计算轴 - 忽略变量定义，因为当成了fft
Nx = length(x); % 信号长度
fs = 1/srs_spacing;
f = (0:Nx-1)*(fs/Nx); % 频率轴

% 绘制谱
figure;
plot(f/TC, abs(X));
grid on;
xlim([0 256])

% 应该是 64 TC，但只能得到60 TC。

%% 复数信号案例-测试时延分辨率
% 创建一个复数信号
% TC = 1/(480 * 1000 * 4096); 

fs = 1/srs_spacing; % 采样频率（Hz）
t = (-num_srs_subcarriers/2:num_srs_subcarriers/2-1)/fs; % 时间向量
f1 = 64*TC; % 信号频率1（Hz）
f2 = 120*TC; % 信号频率2（Hz）
x = exp(1i*2*pi*f1*t) + exp(1i*2*pi*f2*t); % 复数信号

% 对复数信号进行傅里叶变换
X = fft(x);

% 计算频率轴
Nx = length(x); % 信号长度
f = (0:Nx-1)*(fs/Nx); % 频率轴

% 绘制复数信号的频谱
figure;
plot(f/TC, abs(X));
grid on;
xlim([0 256])

% 时延分辨率低为 20TC。

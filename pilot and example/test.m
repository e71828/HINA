%% 相群移
pilot = load("pilot.mat");
example_64Tc = load("example_64Tc.mat");
Xf = pilot.pilot;
Yf = example_64Tc.example_64Tc;

% SRS的频率间隔由下一 section 求得。
TC = 1/(480 * 1000 * 4096); 
tau = -diff(unwrap(angle(Yf./Xf))) / (2*pi*srs_spacing) / TC;  % 计算相群移

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

% 计算频率中心
srs_subcarriers = floor((N - num_srs_subcarriers) / 2) + (1:num_srs_subcarriers);
frequencies = fC + (srs_subcarriers - N/2) * delta_f;

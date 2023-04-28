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
target_length = 50;

grp1 = [126   141   221   238   290];  % 判定过大，大于190。
grp2 = 142;  % 幅值太弱。
%% 处理前400个文件
for i = 1:20
    filename = fullfile(matFiles(i).folder, matFiles(i).name);  % 获取文件名及路径
    data = load(filename);  % 加载MAT文件中的数据
    variable_names = who('-file', filename);  % 获取MAT文件中的变量名
    variable_name = variable_names{1};  % 假设MAT文件中只有一个变量
    Yf = data.(variable_name);  % 获取MAT文件中的变量值    % 对数据进行处理
    Hf = Yf./Xf;
    f_tick_shift = 2*pi*srs_spacing*TC*(1:816);
    % 对原始数据进行移动平均平滑
    smoothed_data = unwrap(angle(Hf));
    for k = 1:10
        smoothed_data = movmean(smoothed_data, 10);
    end
    p = polyfit(f_tick_shift, smoothed_data, 1);

    [~, start_position] = find_decreasing_segment(smoothed_data, target_length);
    p_part = polyfit(f_tick_shift(start_position+(1:target_length)), smoothed_data(start_position+(1:target_length)), 1);

    h = abs(ifft(Hf));
    interp = griddedInterpolant((1:num_srs_subcarriers)*4, h, 'cubic');
    h = interp(1:num_srs_subcarriers*4);

    Nh = length(h); % 信号长度
    h_sel = h(1:52);
    t_sel = (0:51)*(1/srs_spacing/Nh)/TC;

    [result_h, result_t] = find_peaks_with_conditions(h_sel, t_sel);
    
    % 记录第一个延迟峰值
    if isempty(result_t)
        fprintf('No valid peaks found in data of No.%d.\n', i);
    end

    if p(1)/p_part(1) > 5 || p(1)/p_part(1) < 1/5
        if result_t(1) < 50
            tau_ans(i) =  result_t(1);
        else
            tau_ans(i) =  min(-p_part(1),-p(1));
        end
    else
        tau_ans(i) =  result_t(1);
    end



    % 绘制谱
    figure;
    subplot(2,1,1)
    plot(t_sel,h_sel);
    grid on;
    xlim([0 256])
    hold on; stem(result_t, result_h, 'r', 'filled');
    title(['ifft method, delay value: ' num2str(tau_ans(i))])

    subplot(2,1,2)
    plot(f_tick_shift,smoothed_data,'r');
    hold on
    plot(f_tick_shift,polyval(p,f_tick_shift),'g')
    plot(f_tick_shift(start_position+(1:target_length)),polyval(p_part,f_tick_shift(start_position+(1:target_length))), 'b')
    title(['fit with line, p1 = ' num2str(p(1)) ', p_part1 =' num2str(p_part(1))], 'Interpreter','none')
    xlim([0 f_tick_shift(end)])
    plot(f_tick_shift, unwrap(angle(Hf)), 'k')  % 相位复原，难点
    plot(f_tick_shift, 816*abs(Hf), 'm')
    legend('smooth-angle','fit-angle','part-angle','angle','Hf')
    sgtitle(variable_name, 'Interpreter','none');


end 
%% 判断合理
isInRange = (tau_ans(1:400) >= 0) & (tau_ans(1:400) <= 190); % 
allInRange = all(isInRange, 'all');
disp(find(~isInRange));
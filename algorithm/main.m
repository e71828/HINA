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
%% 处理前400个文件
for i = 1:400
    filename = fullfile(matFiles(i).folder, matFiles(i).name);  % 获取文件名及路径
    data = load(filename);  % 加载MAT文件中的数据
    variable_names = who('-file', filename);  % 获取MAT文件中的变量名
    variable_name = variable_names{1};  % 假设MAT文件中只有一个变量
    Yf = data.(variable_name);  % 获取MAT文件中的变量值    % 对数据进行处理
    Hf = Yf./Xf;
    h = abs(ifft(Hf));
    interp = griddedInterpolant((1:num_srs_subcarriers)*4, h, 'cubic');
    h = interp(1:num_srs_subcarriers*4);

    Nh = length(h); % 信号长度
    h_sel = h(1:64);
    t_sel = (0:63)*(1/srs_spacing/Nh)/TC;

    [result_h, result_t] = find_peaks_with_conditions(h_sel, t_sel);
    
    % 记录第一个延迟峰值
    if isempty(result_t)
        fprintf('No valid peaks found in data of No.%d.\n', i);
    else
        tau_ans(i) =  result_t(1);
    end

end 
%% 处理后400个文件
tau_est = zeros(4,1);
for i = 401:800
    filename = fullfile(matFiles(i).folder, matFiles(i).name);  % 获取文件名及路径
    data = load(filename);  % 加载MAT文件中的数据
    variable_names = who('-file', filename);  % 获取MAT文件中的变量名
    variable_name = variable_names{1};  % 假设MAT文件中只有一个变量
    Yf = data.(variable_name);  % 获取MAT文件中的变量值    % 对数据进行处理
    Hf = Yf./Xf;
    f_tick_shift = 2*pi*srs_spacing*TC*(1:816);
    tau_est = zeros(1,4);
    for j=1:size(Hf,1)
        p = polyfit(f_tick_shift, unwrap(angle(Hf(j, :)),[],2), 1);
    
        h = abs(ifft(Hf,[],2));
        interp = griddedInterpolant((1:num_srs_subcarriers)*4, h', 'cubic');
        h = interp(1:num_srs_subcarriers*4)';
    
        Nh = length(h); % 信号长度
        h_sel = h(1:52);
        t_sel = (0:51)*(1/srs_spacing/Nh)/TC;

        [result_h, result_t] = find_peaks_with_conditions(h_sel, t_sel);
        tau_est(j) = result_t(1);
    end
    tau_ans(i) = mean(tau_est);
end 
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
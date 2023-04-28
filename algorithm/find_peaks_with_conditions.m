function [result_h, result_t] = find_peaks_with_conditions(h_sel, t_sel)
    %% 筛选峰值

    % 初始化结果向量
    result_h = [];
    result_t = [];
    
    % 使用findpeaks函数找到所有峰值及其位置
    [peaks, locs] = findpeaks(h_sel);
    % 获取最大峰值
    max_peak = max(peaks);
    % 遍历所有峰值，检查是否满足条件
    for i = 1:length(peaks)
        h = peaks(i);
        t = t_sel(locs(i));
        % 检查是否满足条件
        if (t <= 50 && h > 0.1) || (t > 50 && t <= 100 && h > 0.04) || (t > 100 && h > max_peak/6)
            % 如果满足条件，则将该峰值添加到结果中
            result_h = [result_h; h];
            result_t = [result_t; t];
        end
    end
    if isempty(result_t) && max_peak < 0.04
        result_h = max_peak;
        result_t = t_sel(locs(peaks == max_peak));
    elseif isempty(result_t) && max_peak < 0.06
        result_h = peaks(peaks > max_peak/3);
        result_t = t_sel(locs(peaks > max_peak/3));
    end

end
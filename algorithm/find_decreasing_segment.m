function [decreasing_segment, start_position] = find_decreasing_segment(smoothed_data, target_length)
    % 选择平滑数据的后半段
    half_index = ceil(length(smoothed_data) / 2);
    data_second_half = smoothed_data(half_index:end);

    % 计算一阶差分
    diff_data = diff(data_second_half);

    % 初始化下降段的起始和结束索引
    start_index = 1;
    end_index = 1;

    % 在后半段数据中找出完全下降的一段数据
    for i = 1:length(diff_data)
        if diff_data(i) < 0
            end_index = i + 1;
            % 如果找到的下降段长度达到目标长度，则提前退出循环
            if end_index - start_index + 1 >= target_length
                break;
            end
        else
            start_index = i + 1;
            end_index = start_index;
        end
    end

    % 提取完全下降的一段数据
    decreasing_segment = data_second_half(start_index:end_index);

    % 计算该段在 smoothed_data 中的起始位置
    start_position = half_index + start_index - 1;
end

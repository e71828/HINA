function mergeAnswerFiles(file1, file2, outFile)
% 合并两个答案文件，将结果写入新文件中

if file1 == -1 || file2 == -1
    error('无法打开文件');
end

% 读取并写入第一个文件的内容
while ~feof(file1)
    line = fgetl(file1);
    fprintf(outFile, '%s\n', line);
end

% 读取并写入第二个文件的内容
while ~feof(file2)
    line = fgetl(file2);
    fprintf(outFile, '%s\n', line);
end

fclose(file1);  % 关闭第一个文件
fclose(file2);  % 关闭第二个文件
fclose(outFile);  % 关闭输出文件
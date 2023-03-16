% for s = 1:16
% 
% 
%     textFileName = ['sim' num2str(s) '.txt'];
%     if isfile(textFileName)
%         fid = fopen(textFileName, 'rt');
%         text


clear 
close all

    n_line = 0;
    textFileName = ['sim_test.txt'];
    if isfile(textFileName)
        fid = fopen(textFileName, 'rt');
        line_ex = fgetl(fid);
        while line ~= -1
            line = fgetl(fid);
            n_line = n_line+1;
        end
        frewind(fid);
        sim = zeros(n_line-1,13);
        line_ex = fgetl(fid); %prima linea con i caratteri
        for f=1:n_line-1
            line = fgetl(fid);
            sim(f,:) = str2num(line);
        end
        fclose(fid);
    else
        fprintf('File %s does not exist.\n', textFileName);
    end

    

%     line2 = str2num(line_ex2);
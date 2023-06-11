function Labels = getLabels(path_to_data,fileName)

fid = fopen([path_to_data fileName]);
data = char(fread(fid)');
fclose(fid);

idx = [1,strfind(data,newline)];

Labels = cell(length(idx)-1,1);
for n = 1:length(idx)-1
    if n == 1
        Labels{n} = data(idx(n):idx(n+1)-2);
    else
        Labels{n} = data(idx(n)+1:idx(n+1)-2);
    end
end





function [T,nL,CM] = readTreeFile(path_to_data,fileName)

fid = fopen([path_to_data fileName]);
data = char(fread(fid)');
fclose(fid);

T_start = strfind(data,'(');
T_end = strfind(data,';');
T = data(T_start:T_end-1);

nL = countLeafs(T);
CM = buildTreeData(T,nL);


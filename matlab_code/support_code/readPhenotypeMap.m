function [Y,pi_state,Labels] = readPhenotypeMap(path_to_data,fileName)

fid = fopen([path_to_data fileName]);
data = char(fread(fid)');
fclose(fid);

phenotypeMap = str2num(data(15:end)); %#ok<ST2NM>
[Y,pi_state,Labels] = recodePhenotype(phenotypeMap);
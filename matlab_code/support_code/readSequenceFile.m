function data = readSequenceFile(path_to_data,fileName)

fid = fopen([path_to_data fileName]);
data = char(fread(fid)');
fclose(fid);

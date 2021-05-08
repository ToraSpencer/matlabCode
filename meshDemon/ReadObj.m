%% ��OBJ�ļ��ж�ȡ��������
function [vertexArray] = ReadObj(fullfilename)

fileId = fopen(fullfilename, 'r');
numVertex = 0;
while ~feof(fileId)
    line = fgetl(fileId);
    if line(1) == 'v' 
        numVertex = numVertex + 1; 
        vertexArray(numVertex, :) = GetData(line); 
    end
end
fclose(fileId);
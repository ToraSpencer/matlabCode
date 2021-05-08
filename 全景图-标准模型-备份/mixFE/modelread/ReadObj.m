function [vertexArray] = ReadObj(fullfilename)
% obtain pointArray and faceArray of the given .obj file
%inputs: 
%     fullfilename  -  the full file path
%outputs:
%      vertexArray  -  point array
%        faceArray  -  face array

fileId = fopen(fullfilename, 'r');
numVertex = 0;
numFace = 0;
while ~feof(fileId)
    line = fgetl(fileId);
    if line(1) == 'v' 
        numVertex = numVertex + 1; 
        vertexArray(numVertex, :) = GetData(line); 
    end
    if line(1) == 'f' 
        numFace = numFace + 1;
        faceArray(numFace, :) = GetData(line); 
    end 
end
fclose(fileId);
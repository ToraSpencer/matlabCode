function [vertex,face] = ReadObj(fullfilename)
% obtain pointArray and faceArray of the given .obj file
%inputs: 
%     fullfilename  -  filename
%outputs:
%      vertexArray  -  point array
%        faceArray  -  face array
%
%---------------------------------------------------------
% fullfilename = 'output_roots.obj';
fileId = fopen(fullfilename, 'r');
numVertex = 0;
numFace = 0;
num = 1;
Line = char();
N = 1;
while ~feof(fileId)
    
    line = fgetl(fileId);
   
    if line(1) == 'v' 
        Line = [Line;line(1:2)];
        num = num+1;
        numVertex = numVertex + 1; 
        vertexArray(numVertex, :) = GetData(line); 
        
    end
    if line(1) == 'f' 
        Line = [Line;line(1:2)];
        
            
        numFace = numFace + 1;
        faceArray(numFace, :) = GetData(line);
        num = num+1;
    end
    
    if feof(fileId)
        vertex = vertexArray;face = faceArray;
    end
           
end
fclose(fileId);
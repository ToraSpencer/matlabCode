%% 从OBJ文件中读取网格数据

function [all_struct] = Read_Obj(fullfilename)
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
    if num >2
        if Line(num-1,1) == 'v'&&Line(num-2,1) == 'f'
            all_struct(N,1).vertex = vertexArray(1:end-1,:);all_struct(N,1).face = faceArray;
            %reset
            numVertex = 1;
            numFace = 0;
            vex = vertexArray(end,:);
            vertexArray = vex;
            faceArray = [];
            N = N+1;
        end
    end
    if feof(fileId)
        all_struct(N,1).vertex = vertexArray;all_struct(N,1).face = faceArray;
    end
           
end
fclose(fileId);
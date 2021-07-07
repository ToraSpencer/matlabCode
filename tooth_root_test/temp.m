%% 提取测试数据

clc;
clear all;

load('D:\workstation\gitRepositories\matlabCode\tooth_root_test\matlab.mat')
toothdata = handle.model.toothdata;
nT = 14;

toothCenter = zeros(nT,3);
for i = 1:nT
    toothCenter(i,:) = mean(toothdata{1,i});
end
dentalCenter = mean(toothCenter);

fdi = 0;
for i = 1:nT

    if(i < 8)
        fdi = 48 - i;
    else
        fdi = 23 + i;
    end
    
    
    V0 = toothdata{1,i};
    F0 = toothdata{2,i};

    % 提取待拼接牙齿的牙冠区域
    patientCutVers = V0(1:toothdata{8,i},:);
    patientCutTris = F0(sum(F0 <= toothdata{8,i}, 2)  == 3, :);
    
    nameStr = ['patientTooth', num2str(fdi), '.obj'];
    writeOBJ(nameStr, patientCutVers, patientCutTris);
    
    patientDir = toothdata{7,i};
    patientAxisTrans(:,3) = -patientDir';
    patientAxisTrans(:,1) = normalizerow(cross(toothCenter(i,:) - dentalCenter, -patientDir))';
    patientAxisTrans(:,2) = cross(patientAxisTrans(:,3)', patientAxisTrans(:,1));
    
    nameStr = ['axisPatient', num2str(fdi), '.obj'];
    OBJwriteVertices(nameStr, patientAxisTrans');
    
 
    n = 20;

    hole = Calc_Boundary(patientCutTris);
    gumlineIdx = hole.boundary.edge(:,1);
    gumline = patientCutVers(gumlineIdx,:);
    nameStr = ['gumline', num2str(fdi), '.obj'];
    OBJwriteVertices(nameStr, gumline);
 
end
        disp('finished');

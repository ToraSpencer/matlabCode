clc;
clear all;


%%
rEdgeVersIdx = [3383
3699
4083
3938
4308
4099
3860
3767
3393
4004
4109
3975
3889
3734
4325
4286
3411
4060
3575
3641
3599
3410
3968
3667
3673
3936
4174
3675
4057
3739
4169
3584
4100
4064
3736
3413
3872
3836
4055
3680
4093
3922
4227
4192
3945
3629
4152
4317
3994
4248
4291
3826
3786
4061
4158
];

    [patientCutVers, patientCutTris] = readOBJ('./step5/切割牙冠cpp.obj');
    [finalVers, finalTris] = readOBJ('./step5/融合网格cpp.obj');

   adjMat = adjacency_matrix(finalTris);   % 融合网格的邻接矩阵,维度为versCount*versCount，两个顶点有边连接则对应元素为1；
 
    temp = sum(adjMat, 2);   
    L = adjMat - diag(sparse(sum(adjMat,2)));
    
    %for debug
    flag = (temp == zeros(size(temp, 2), 1));
    
    A = L;
    
    L2 = L*L;
    A(rEdgeVersIdx, :) = L2(rEdgeVersIdx, :);
 
    B = L*finalVers;
    B(rEdgeVersIdx, :) = 0;
    
    Acon = 1:size(patientCutVers,1);
    Bcon = patientCutVers;
    finalVers = solve_equation_modified(A, B, Acon, Bcon);
    
    
    writeOBJ('最终结果.obj', finalVers, finalTris);

    disp('finished');

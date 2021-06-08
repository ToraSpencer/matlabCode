function DatWriteIdxs(filename, matIdx)
f = fopen( filename, 'w' );

% 前两行是行列信息：
fprintf(f, 'row %d\n', size(matIdx, 1));
fprintf(f, 'col %d\n', size(matIdx, 2));

tempMat = matIdx - ones(size(matIdx));          % 写到文件里，索引改为从0开始，而不是matlab里从1开始。
fprintf( f, '%d\n', tempMat);

 

fclose(f);
function DatWriteIdxs(filename, matIdx)
f = fopen( filename, 'w' );

% ǰ������������Ϣ��
fprintf(f, 'row %d\n', size(matIdx, 1));
fprintf(f, 'col %d\n', size(matIdx, 2));

tempMat = matIdx - ones(size(matIdx));          % д���ļ��������Ϊ��0��ʼ��������matlab���1��ʼ��
fprintf( f, '%d\n', tempMat);

 

fclose(f);
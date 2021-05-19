function vecWriteUnsigned( fileName, vec )

m1 = ones(size(vec));
f = fopen( fileName, 'w' );

fprintf( f, '%d\n', (vec - m1)');

 

fclose(f);



end


function OBJwriteVertices(filename, Vers)
f = fopen( filename, 'w' );

fprintf( f, 'v %0.17g %0.17g %0.17g\n', Vers');

 

fclose(f);

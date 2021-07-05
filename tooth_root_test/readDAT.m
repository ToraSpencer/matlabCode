function [ indexes] = readDAT( filename )
 
  V = zeros(10000,3);
  F = zeros(10000,ss);
  UV = zeros(10000,3);
  TF = zeros(10000,ss);
  N = zeros(10000,3);
  NF = zeros(10000,ss);

  triangulated = false;
  all_ss = true;
  fp = fopen( filename, 'r' );
  type = fscanf( fp, '%s', 1 );
  count = 0;
  while strcmp( type, '' ) == 0
      line = fgets(fp);
      if strcmp( type, 'v' ) == 1
          v = sscanf( line, '%lf %lf %lf' );
          numv = numv+1;
          if(numv>size(V,1))
            V = cat(1,V,zeros(10000,3));
          end
          V(numv,:) = [v(1:3)'];
      elseif strcmp( type, 'vt')
          v = sscanf( line, '%f %f %f' );
          numuv = numuv+1;
          if size(UV,2)>2 && length(v) == 2
              UV = UV(:,1:2);
          end
          if(numuv>size(UV,1))
            UV = cat(1,UV,zeros(10000,length(v)));
          end
          UV(numuv,:) = [v'];
      elseif strcmp( type, 'vn')
          n = sscanf( line, '%f %f %f' );
          numn = numn+1;
          if(numn>size(N,1))
            N = cat(1,N,zeros(10000,3));
          end
          N(numn,:) = [n'];
      elseif strcmp( type, 'f' ) == 1
          [t, count] = sscanf(line,'%d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d');
          if (count>2)
              tf = t(2:3:end);
              nf = t(3:3:end);
              t = t(1:3:end);
          else
              [t, count] = sscanf(line, '%d/%d %d/%d %d/%d %d/%d %d/%d');
              if (count>2)
                  tf = t(2:2:end);
                  t = t(1:2:end);
                  nf = -ones(numel(t),1);
              else
                [t, count] = sscanf(line, '%d//%d %d//%d %d//%d %d//%d %d//%d');
                if (count>2)
                    nf = t(2:2:end);
                    t = t(1:2:end);
                    tf = -ones(numel(t),1);
                else
                    [t, count] = sscanf( line, '%d %d %d %d %d %d %d %d %d %d %d\n' );
                    tf = -ones(numel(t),1);
                    nf = -ones(numel(t),1);
                end
              end
          end
          assert(numel(t) == numel(tf));
          assert(numel(t) == numel(nf));
          if numel(t) > ss
            if ~triangulated
              warning('Trivially triangulating high degree facets');
            end
            triangulated = true;
          end
          j = 2;
          i = 1;
 
          while true
            if numel(t) > ss
              corners = [1 2 3];
 
            else
              if all_ss && numel(t)<ss
                warning('Small degree facet found');
                all_ss = false;
              end
              corners = 1:numel(t);
            end
            numf = numf+1;
            if(numf>size(F,1))
              F = cat(1,F,zeros(10000,ss));
            end
            F(numf,1:numel(corners)) = [t(corners)'];
            numtf = numtf+1;
            if(numtf>size(TF,1))
              TF = cat(1,TF,zeros(10000,ss));
            end
            TF(numtf,1:numel(corners)) = [tf(corners)'];
            numnf = numnf+1;
            if(numnf>size(NF,1))
              NF = cat(1,NF,zeros(10000,ss));
            end
            NF(numnf,1:numel(corners)) = [nf(corners)'];
            if numel(t) <= ss
              break;
            end
            t = t([1 3:end]);
            %VtA = VtA([1 3:end],:);
            if numel(t) < 3
              break;
            end
          end
      elseif strcmp( type, '#' ) == 1
      end
      type = fscanf( fp, '%s', 1 );
  end
  fclose( fp );
 
  V = V(1:numv,:);
  F = F(1:numf,:);
  UV = UV(1:numuv,:);
  TF = TF(1:numtf,:);
  N = N(1:numn,:);
  NF = NF(1:numnf,:);
 
  if (size(UV,1)>0) UV = UV; end

end

 


 


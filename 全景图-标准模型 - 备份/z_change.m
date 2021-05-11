function V_changed  = z_change(V, y,thcenter,p_apx,rat,va)
l_y = (y(floor(length(y)/2)+1) - mean([y(1),y(end)]))*rat;
l =  thcenter(3) - p_apx(3);
p_apx_changed(3) = l_y+thcenter(3);
if va == 1     %varargin= 1:ÏÂò¢
    n = find (V(:,3)<thcenter(3));
else
    n = find (V(:,3)>thcenter(3));
end
V(n,3) = (V(n,3) - thcenter(3))*l_y/l+thcenter(3);
V_changed = V;
     


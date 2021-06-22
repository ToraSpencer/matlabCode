function V = smooth_loop(V, w)
    nV = size(V,1);
    
    vl = mod((1:nV)-2, nV) + 1;
    vr = mod(1:nV, nV) + 1;
    W = sparse([1:nV, 1:nV], [vl, vr], 1);
    A1 = speye(nV) - spdiags(1./sum(W,2), 0, nV, nV) * W;
	B1 = zeros(size(V));
    A2 = speye(nV);
    B2 = V;
	
    V = [A1; w*A2]\[B1; w*B2];
end

function [ versOut ] = mySmooth( versIn, param )
    versCount = size(versIn,1);
    
    vl = mod((1:versCount)-2, versCount) + 1;                 % mod()ศกำเสý
    vr = mod(1:versCount, versCount) + 1;
    W = sparse([1:versCount, 1:versCount], [vl, vr], 1);
    A1 = speye(versCount) - spdiags(1./sum(W,2), 0, versCount, versCount) * W;
	B1 = zeros(size(versIn));
    A2 = speye(versCount);
    B2 = versIn;
	
    versOut = [A1; param*A2]\[B1; param*B2];
end
 

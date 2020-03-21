function varargout = myLegendre_A2(A1,np,P)
% Start from 2
% [cAA Yf] or cAA = myLegendre_A(Y, n, P);
% Yf = (cAA(:)'* P')';
% P[DIM1 DIM2], DIM1 = points, DIM2 = number of facts
nl = size(P, 2); cAA = zeros(nl, 1);
A0 = (A1(1:(end-1))+A1(2:end))/2;
if np>nl; np=nl; end;
nf = 1./(size(P,1)-1);    
for ii=1:np
    PP = P(:,ii);
    %PP = (PP(1:(end-1))+PP(2:end))/2; %PP = mean([PP(1:(end-1)) PP(2:end)], 2);
    cAA(ii)=sum(PP(:).*A0(:)*(2*ii-1)) * nf;
end
varargout{1}=cAA(1:np);
if nargout > 1
    varargout{2}=(cAA(:)'* P')';
end
end
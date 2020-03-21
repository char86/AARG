function [Ao Eo]=myGauss1D_fit2(X,Y, At, Ft)
% Using: [A E] = myGauss1D_fit2(X,Y); same as myGauss1D_fit, used for Negtive X fitting;
% A =[X0 S Offset Amp];
Ao=[0 0 0 0]';Eo=[0 0 0 0]'; A=[0 0 0 0]';E=[0 0 0 0]';
mm=[max(Y(:)) min(X(:))]; if mm(1)==0; return; else Y=Y/mm(1); end; X=X-mm(2);
A(3) = min(Y(:));
A(4) = (1 - A(3))*1.1;
if abs(sum(Y - A(3)))<1e-7; Y=rand(size(Y))*1e-7;end;
A(1) = sum((Y - A(3)).* X)/sum(Y - A(3));
nIndex = find(Y>mean(Y(:)));%nIndex = find(Y>A(1));
if numel(nIndex)>0
   A(2) = (max(X(nIndex)) - min(X(nIndex)))/2;
else
    A(2) = 2;
end
F=[1 1 1 1]';
if nargin>3;
    ii=find(Ft==0); F(ii)=0; A(ii)=At(ii);
end
[L C]=myGAUSS_FIT(X, Y, A, F, E);
A(4)=A(4)*mm(1); A(1)=A(1)+mm(2); Ao = A; Eo = E;
end

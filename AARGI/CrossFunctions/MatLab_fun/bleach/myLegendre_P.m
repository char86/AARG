function [P X]=myLegendre_P(varargin)
% Legendre [-1 1]
% Usage: [P X]=myLegendre_P(nl);
% nl=[1000 1000]; % points,  number of facts
if nargin==0; nl=[1000 1000];
else nl=varargin{1};end
if(nl(2)<5);nl(2)=5;end
X=(-1:(2/(nl(1)-1)):1)';X(1)=-1;X(end)=1;
P=zeros(nl);
P(:,1)=X.*0+1;
P(:,2)=X;
for ii=3:nl(2)
    n=ii-2;
    P(:,ii)=((2*n+1).*X.*P(:,ii-1)-n*P(:,ii-2))/(n+1);
end
end

function [Out]=myCros_00_d3(A, B, nDim)
% Using:
%   myCros_00_d3(aIN_Large, aIN_Small[, dim]); aIN is input data, dim, by default, is the last one.
if nargin<3; nDim=[]; end;
ss=size(A); if numel(nDim)<1||nDim>numel(ss); nDim=numel(ss);end
% following 2 line is OK, but last part is zeros
%if nDim~=1; ii=1:numel(ss); ii(nDim)=1; ii(1)=nDim; A=permute(A,ii); end;
% Out=myCros_00_d3_vf(A,B,1);

if nDim~=1; ii=ss(nDim); ss(nDim)=ss(1); ss(1)=ii; ii=1:numel(ss); ii(nDim)=1; ii(1)=nDim;A=permute(A,ii); end;
Out=myCros_00(A,B); Out=reshape(Out,ss);

if nDim~=1; Out=permute(Out,ii); end;
end
function [varargout]= myReSize3_bin(dataS, nn, bNorm)
% function [Out s2]= myReSize3_bin(dataS, nn, bNorm)
% [Img Size]= myReSize3_bin(Img_3d, bin_size)
% resize IMG, repeat in 3ed DIM
% cf myReSize3, using fft
if nargin<2; nn=3; elseif numel(nn)<1; nn=3; end
if nn>0
    nn=floor(nn); if nn<1; nn=1; end
    if nn==1; Out=dataS;
    else % for 3d, cut and resize
        ss=size(dataS); dataS2 = dataS(1:(floor(ss(1)/nn)*nn),1:(floor(ss(2)/nn)*nn),:);
        s2=mySize(dataS2); dataS2=squeeze(sum(reshape(dataS2,[s2(1),nn, s2(2)/nn,s2(3)]),2));
        s2=mySize(dataS2); dataS2=squeeze(sum(reshape(dataS2,[nn,s2(1)/nn,s2(2),s2(3)]),1));
        Out=dataS2;
    end
else % expend
    nn=floor(0-nn); if nargin<3; bNorm=0; end;
    if nn<=1; Out=dataS;
    else % expend in 3D
        ss=mySize(dataS); Out=zeros([ss(1)*ss(2)*nn*nn ss(3)]);
        dataS=reshape(dataS, [], ss(3));
        xx = (1:ss(1))-1; xx=xx'*nn;
        xx = repmat(xx, [1 ss(2)])+repmat(((1:ss(2))-1)*ss(1)*nn*nn, [ss(1) 1])+1;
        xx = reshape(xx, [], 1);
        for ii=1:nn; for jj=1:nn;
                kk=(ii-1)+(jj-1)*ss(1)*nn;
                Out(xx+kk,:)=dataS;
            end; end
        Out=reshape(Out,[[ss(1) ss(2)]*nn ss(3)]);
    end
end
varargout{1}=Out; if nargout>1;s2=size(dataS2); varargout{2}=s2; end
    function ex = mySize(In)
        ex=size(In); if numel(ex)<3;ex(3)=1;end
    end
end

% Following is test code for -bin
% B=myReSize3_bin(A, -2);
% %%
% A=reshape(1:200, [10 20]); B=zeros(20,40);
% ss=[10,20]; nn= 2;
% 
% xx = (1:ss(1))-1; xx=xx'*nn;
% xx = repmat(xx, [1 ss(2)])+repmat(((1:ss(2))-1)*ss(1)*nn*nn, [ss(1) 1])+1;
% mm=zeros(nn);
% for ii=1:nn; for jj=1:nn;
%         kk=(ii-1)+(jj-1)*ss(1)*nn;
%         B(xx+kk)=A(:);
%     end; end

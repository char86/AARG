function varargout = myNorm_G(dd, AA, opt)
% Using:
% Out = myNorm_G(data, AA, [dim, mean_A(1), nofig, mean]);
% data default is last dim
% AA shall be the same as data except the data-dim
% opt(2)=0; mean value will be from opt(4)
%
%
% e.g.
% AA = myGauss1D_All(data);
% data = myNorm_G(data, AA, [-1 1 1]); % same as data = myNorm_G(data, AA);
%
% default dim=-1, last dim
% default mean_A1, take all AA amp as global amp
%
%
ss0=size(dd); nn=numel(ss0);
if nargin<3; opt=[]; end; if numel(opt)<1; opt=[-1 1 1]; end; if numel(opt)<6; opt(6)=0; end
if opt(1)<1||opt(1)>nn; dim=nn; else dim=opt(1); end

if dim~=nn; pp=1:nn; pp(dim)=nn; pp(nn)=dim; dd=permute(dd,pp); AA=permute(AA,pp); end;
dd=reshape(dd,[],ss0(dim)); AA=reshape(AA,[],4); smm=mean(AA);  smm=smm(1); mm=size(dd,1); Ao=zeros(mm,ss0(dim)); %Ao=zeros(ss0(dim),mm);
if opt(2)==0; smm=opt(4); end
if opt(3)==0; set(gcbf,'pointer','watch'); hh = waitbar(0,'Normalizing ... '); set(hh,'Name','Waiting for Normalization.'); llh = ceil(mm/10);end
for ii=1:mm;
    kk=double(squeeze(dd(ii,:))); A=AA(ii,:);
    if A(2)>1e-3; kk=(kk-A(1))/A(2)+smm; else kk= kk-A(1) + smm; end; Ao(ii, :)=kk; %Ao(:,ii)=kk;
    if opt(3)==0; if mod (ii,llh)==1;  waitbar(floor(ii/llh)/10,hh, ['Normalizing... ' num2str(ii) '/' num2str(mm)]); end; end
end;
if opt(3)==0; if ishandle(hh); close(hh); end; set(gcbf,'pointer','arrow');end
%Ao=Ao'; 
Ao=reshape(Ao, ss0);
if dim~=nn; Ao=permute(Ao,pp); end;
varargout{1}=Ao;
end

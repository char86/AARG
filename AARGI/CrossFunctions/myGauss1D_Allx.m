function [Ao]=myGauss1D_Allx(dd, dim, ll, noFig)
% Using: [Ao]=myGauss1D_Allx(dd, dim, ll, noFig); % Ao =[X0 S Offset Amp];
%
ss=size(dd); nn=numel(ss);
if nargin<2;dim=[];end; if numel(dim)<1; dim=nn; end
if nargin<3;ll=[];end; if numel(ll)<1;ll=50;end
if nargin > 3;  bnoFig=1; else  bnoFig=0;  noFig=0; end

if dim~=nn; pp=1:nn; pp(dim)=nn; pp(nn)=dim; dd=permute(dd,pp); end;
ss0=size(dd); dd=reshape(dd,[],ss0(nn)); mm=size(dd,1); Ao=zeros(mm,4); tt=[-ll:ll];
if bnoFig==0
    set(gcbf,'pointer','watch');
    hh = waitbar(0,'Distribution Gauss Fitting ... '); set(hh,'Name','Waiting for Distribution Gauss Fitting.'); llh = ceil(mm/10);
end
%tt=[min(dd(:)) max(dd(:))]; tt(3)=(tt(2)-tt(1))./4096; ttl=tt(1):tt(3):tt(2);
%if tt(3)>1e-7;
    for ii=1:mm;
        kk=squeeze(dd(ii,:)); ttl=[min(kk(:)) max(kk(:))]; ttl(3)=(ttl(2)-ttl(1))./ll/2.5; ttl=ttl(1):ttl(3):ttl(2);
        [Y X]=hist(kk,ttl); [jj jm]=max(Y); if (X(jm)-2049)<1&&jm>1&&jm<numel(X); Y(jm)= (Y(jm-1)+Y(jm+1))/2; end; %Y(2049)=(Y(2048)+Y(2050))/2;
        [jj jm]=max(Y); jj=jm+[-ll:ll]; if jj(1)<1;jj=jj-jj(1)+1;end;
        if jj(ll*2+1)>numel(Y); jj=jj-jj(ll*2+1)+numel(Y); jm=find(jj>1); jj=jj(jm); end
        A=myGauss1D_fit2(X(jj),Y(jj)); Ao(ii,:)=A;
        if noFig==2; figure;plot(X,Y);  hold all; plot(X, myGauss1D_Cal(X,A));xlim([A(1)-17*A(2),A(1)+17*A(2)]); end %A =[X0 S Offset Amp];
        if bnoFig==0
            if mod (ii,llh)==1;  waitbar(floor(ii/llh)/10,hh, ['Distribution Gauss Fitting ... ' num2str(ii) '/' num2str(mm)]);       end;
        end
    end;
%end
if bnoFig==0; if ishandle(hh); close(hh); end; set(gcbf,'pointer','arrow');end
Ao=reshape(Ao,[ss0(1:(nn-1)) 4]);
if dim~=nn; Ao=permute(Ao,pp); end;

end

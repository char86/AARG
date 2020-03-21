function Out = myEXP_Curves(ll, tt, ttn)
% Using:
% Out = myEXP_Curves([Zero_Before EXP_Length], [List_of_Tau], [List_of_Tau_Neg])
% Test Code:
% ll=[7 17];tt=[250 200 170 140 110 75 50 30 20 15 10 7 5 3 2 1];
% tmp = myEXP_Curves(ll, tt);
% Test Code2:
% ll=[7 17];tt=[10]; ttn=[0.1 0.3 0.5 0.7 1 1.3 1.7 2.2 3];
% tmp = myEXP_Curves(ll, tt, ttn);
xx=1:ll(2);
if nargin<3
    ee=zeros(numel(tt),ll(1)+ll(2)); for ii=1:numel(tt);yy=exp(-(xx-1)/tt(ii));ee(ii,ll(1)+(1:ll(2)))=yy;end
else
    if numel(tt)==numel(ttn)
        ee=zeros(numel(tt),ll(1)+ll(2)); for ii=1:numel(tt); if tt(ii)<=ttn(ii); continue; end; yy=exp(-(xx-1)/tt(ii)) - exp(-(xx-1)/ttn(ii));ee(ii,ll(1)+(1:ll(2)))=yy;end
    else
        ee=zeros(numel(ttn),ll(1)+ll(2)); for ii=1:numel(ttn);if tt(1)<=ttn(ii); continue; end; yy=exp(-(xx-1)/tt(1)) - exp(-(xx-1)/ttn(ii));ee(ii,ll(1)+(1:ll(2)))=yy;end
    end
end
Out=ee;
end
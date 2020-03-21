function [mm mm0] = myNB3(dd, opt)
% Using:
% Out=myNB3(In);
% or [Out Out0]=myNB3(In, [thread_number, norm]);
% 


if nargin<2; opt=[]; end; if numel(opt)<1;opt=[3 1];end; opt(10)=0;

ss=size(dd); nn=numel(ss);

if nn<3; mm=zeros(ss(1:2)); mm0=mm; return; end
if nn>=4; mm0=nb3(dd, opt(1)); else 
    mm0=nb3(reshape(dd,[ss(1) ss(2) 1 ss(3)]),opt(1));
end
if opt(2)==1; mm=mm0-min(mm0(:));mm=mm./max(mm(:)); else mm=mm0; end
end
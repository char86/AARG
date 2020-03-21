function varargout = myBleachingL3(varargin)
% Using:
%    [Out Bleach]= myBleachingL3(In, n);
% or [Out] = myBleachingL3(In);
%
% Here:
% In: 3-D array (e.g. 2D, time series);
% n:  2 ~ 7, factors needs to remove, by default, n = 3
%
% On return:
% Out: output array with bleaching part removed
% Bleach: the bleaching part
%
% Note: Out more or less ~= In ./ Bleach .* repmat( Bleach(:,:,1), [1 1 size(In, 3)]);
%
% The function removes a single decay part of a fluorescence based on Legendre polynomials.
% For detail info, please refer: Bao & Schild (2014) Fast and Accurate Fitting and Filtering
% of Noisy Exponentials in Legendre Space. PLoS ONE 9(3): e90500. doi:10.1371/journal.pone.0090500
%
% Any question or bug report, please contact us, Thanks. gbao@gwdg.de
%

if nargin > 1;
    nn = varargin{2};if nn<2; nn=2; end
    if nn>17; nn = 17; end
else nn = 3; end
In = varargin{1}; Out=[]; ll=size(In);
if numel(ll)~=3;  fprintf('Please use the function for 3-D array (e.g. time serious 2D data)\n'); return; end
% Start
Out = reshape(In, [], ll(3)); kk=ll(1)*ll(2);
set(gcbf,'pointer','watch');
hh = waitbar(0,'Bleach Correcting ... '); set(hh,'Name','Waiting for Image Bleach Correction.');
% addpath
functionname='myBleachingL3.m';
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath(functiondir);
addpath([functiondir 'bleach']);
% get P
nl=[ll(3) nn];[P X]=myLegendre_P(nl); P=P(:,1:nn);P2=P'; P = (P(1:(end-1),:)+P(2:end,:))/2;
Off = 0; Out = Out - Off; mm=min(Out(:));
if mm<=0;    Off = Off + mm; Out = Out - mm + 1e-7; end
mm = ceil(kk / 10);
if nargout > 1
    Out2 = Out;
    for ii=1:kk
        Y = Out(ii,:); cAA = myLegendre_A2(Y, nn, P)';Yf = cAA* P2;
        Out(ii,:) = Y(:)./Yf(:).*Yf(1); Out2(ii,:) = Yf(:);
        if mod ( ii, mm )==1;  waitbar(floor(ii/mm)/10,hh, ['Bleach Correcting ... ' num2str(floor(ii/kk*100)) '%']);       end;
    end
else
    for ii=1:kk
        Y = Out(ii,:); cAA = myLegendre_A2(Y, nn, P)';Yf = cAA* P2;
        Out(ii,:) = Y(:)./Yf(:).*Yf(1);
        if mod ( ii, mm )==1;  waitbar(floor(ii/mm)/10,hh, ['Bleach Correcting ... ' num2str(floor(ii/kk*100)) '%']);       end;
    end
end
if ishandle(hh); close(hh); end
set(gcbf,'pointer','arrow');
Out=reshape(Out, ll) + Off; mm=[min(In(:)) max(In(:))]; mm(3) = (mm(2)-mm(1))*0.1; mm(1)=mm(1)-mm(3); mm(2)=mm(2)+mm(3); Out(Out<mm(1))=mm(1); Out(Out>mm(2))=mm(2);
varargout{1} = Out; %varargout{1} = reshape(Out, ll) + Off;
if nargout > 1
    varargout{2}=reshape(Out2, ll);
end
end

function Out = myShiftCorrect2D(varargin)
% Using:
%    Out = myShiftCorrect2D(In, s);         % fast method
%    Out = myShiftCorrect2D(In, s, 'c');    % cycle shift
%    Out = myShiftCorrect2D(In, s, 'fft');  % fft shift
if nargin > 1
    data = varargin{1};s = varargin{2};
else
    fprintf('using: Out = myShiftCorrect2D(In, s [, ''c''/''fft'']); \n');return;
end
if nargin > 2;   mm = myName(varargin{3});  else mm = 0; end
if nargin > 3;  bFig=1; else  bFig=0;  end

set(gcbf,'pointer','watch');
hh = waitbar(0,'Shift Correcting ... '); set(hh,'Name','Waiting for Image Shift Correction.');


OneS=zeros(size(data)); jj = size(data,3); ll = ceil(jj/10);
if mm==0; % Correct the Shift
    for ii=1:jj;OneS(:,:,ii)=myShift1_vF(double(data(:,:,ii)), [s(ii,2) s(ii,1)]);
        if mod ( ii, ll )==1;  waitbar(floor(ii/ll)/10,hh, ['Shift Correcting ... ' num2str(ii) '/' num2str(jj)]);       end;
    end;
elseif mm==2; % Correct the Shift with fft
    for ii=1:jj;OneS(:,:,ii)=myShift1_fft(double(data(:,:,ii)), [s(ii,2) s(ii,1)]);
        if mod ( ii, ll )==1;  waitbar(floor(ii/ll)/10,hh, ['Shift Correcting ... ' num2str(ii) '/' num2str(jj)]);       end;
%        if mod ( ii, ll )==1; fprintf('%d ', 10 - floor(ii/ll));end;
    end;
elseif mm==1; % Correct the Shift with cycle shift
    for ii=1:jj;OneS(:,:,ii)=myShift1_vF_CYC(double(data(:,:,ii)), [s(ii,2) s(ii,1)]);
        if mod ( ii, ll )==1;  waitbar(floor(ii/ll)/10,hh, ['Shift Correcting ... ' num2str(ii) '/' num2str(jj)]);       end;
%        if mod ( ii, ll )==1; fprintf('%d ', 10 - floor(ii/ll));end;
    end;
else
    OneS = data;
end
Out = OneS;
if ishandle(hh); close(hh); end
set(gcbf,'pointer','arrow');


    function mmOut = myName(mmIn)
        if mmIn == 'c' ; mmOut = 1;
        elseif mmIn == 'fft'; mmOut = 2;
        else mmOut = 0;end;
    end
end
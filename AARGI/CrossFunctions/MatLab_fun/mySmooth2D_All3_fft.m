
function [Out]=mySmooth2D_All3_fft(data, x, y, wbstr)
% Out = mySmooth2D_All3_fft(In, x, y)
% e.g Out = mySmooth2D_All3_fft(data, 0.7, 0.7);
Out=zeros(size(data));
mn = size(data);ss = [y/2 x/2];
ss(1)=ss(2)*(mn(2)/mn(1)).^2;
[XX YY]=meshgrid((1:mn(2))-mn(2)/2,(1:mn(1))-mn(1)/2);
fGs=exp(-2*pi*pi*((XX*ss(2)/mn(1)).^2+(YY*ss(1)/mn(2)).^2));
set(gcbf,'pointer','watch');
if nargin < 4; wbstr = 1; end
if wbstr == 1
    hh = waitbar(0,'Smoothing ... '); set(hh,'Name','Waiting for Image Smoothing.');
end
jj = size(data,3); ll = ceil(jj/10);
for ii=1:jj;
    Out(:,:,ii)  = abs(ifft2(fftshift(fft2(data(:,:,ii))).*fGs));
    if mod ( ii, ll )==1 && wbstr == 1;  
        waitbar(floor(ii/ll)/10,hh, ['Smoothing ... ' num2str(ii) '/' num2str(jj)]);       
    end;
end
if wbstr == 1; if ishandle(hh); close(hh); end; end
set(gcbf,'pointer','arrow');
end
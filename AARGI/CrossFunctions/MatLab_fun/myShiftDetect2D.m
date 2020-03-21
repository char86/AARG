function Out = myShiftDetect2D(varargin)
% Using:
% Out = myShiftDetect2D(In); % first 1 will be Ref
% Out = myShiftDetect2D(In, Ref);
set(gcbf,'pointer','watch');
hh = waitbar(0,'Shift Detecting ... '); set(hh,'Name','Waiting for Image Shift Detect.');
if nargin == 1
    In = double(varargin{1});
    Out = zeros(size(In, 3), 2); jj=size(In, 3); ll = ceil(jj/10);
    tif_Key = conj(fft2(myPad_2D_in(myFFTNormal_h_in( In(:,:,1) ))));
    Out(1,1) = 0; Out(1,2)=0;
    nHalf = [size(In, 1) + 1, size(In, 2) + 1];
    ss=-6:6;[XI,YI] = meshgrid(ss+7,ss+7);F = [1 1 1 1 1 1]; E = [0 0 0 0 0 0];
    for ii=2:jj
        tifShift = abs(fftshift(ifft2(fft2(myPad_2D_in(myFFTNormal_h_in(In(:,:,ii)))).* tif_Key )));
        [r c] = find(tifShift == max(tifShift(:)));
        tifShiftP=tifShift(ss+r(1),ss+c(1));
        A = [max(tifShiftP(:))-min(tifShiftP(:)),7,1,7,1,min(tifShiftP(:))];
        [Chi n] = Gauss2D_vF2( [XI,YI] , tifShiftP, A, F, E);
        Out(ii,1) = -(A(2)-7+(c-nHalf(2))); % for x
        Out(ii,2) = -(A(4)-7+(r-nHalf(1))); %for y
        if mod ( ii, ll )==1;  waitbar(floor(ii/ll)/10,hh, ['Shift Detecting ... ' num2str(ii) '/' num2str(jj)]);       end;
    end
else
    In = double(varargin{1});
    Out = zeros(size(In, 3), 2); jj=size(In, 3); ll = ceil(jj/10);
    tif_Key = conj(fft2(myPad_2D_in(myFFTNormal_h_in(double(varargin{2})))));
    %Out(1,1) = 0; Out(1,2)=0;
    nHalf = [size(In, 1) + 1, size(In, 2) + 1];
    ss=-6:6;[XI,YI] = meshgrid(ss+7,ss+7);F = [1 1 1 1 1 1]; E = [0 0 0 0 0 0];
    for ii=1:jj
        tifShift = abs(fftshift(ifft2(fft2(myPad_2D_in (myFFTNormal_h_in(In(:,:,ii)))).* tif_Key )));
        [r c] = find(tifShift == max(tifShift(:)));
        tifShiftP=tifShift(ss+r(1),ss+c(1));
        A = [max(tifShiftP(:))-min(tifShiftP(:)),7,1,7,1,min(tifShiftP(:))];
        [Chi n] = Gauss2D_vF2( [XI,YI] , tifShiftP, A, F, E);
        Out(ii,1) = -(A(2)-7+(c-nHalf(2))); % for x
        Out(ii,2) = -(A(4)-7+(r-nHalf(1))); %for y
        if mod ( ii, ll )==1;  waitbar(floor(ii/ll)/10,hh, ['Shift Detecting ... ' num2str(ii) '/' num2str(jj)]);       end;
    end
end
if ishandle(hh); close(hh); end
set(gcbf,'pointer','arrow');
%tifShift = mycrossCorr2D(tif1,tif2);
    function [Out] = myPad_2D_in(In)
        S = size(In);
        Out = zeros(S(1:2)*2);
        Out(1:S(1),1:S(2)) = In(:,:,1);
    end
    function [G]=myFFTNormal_h_in(Img)
        G = zeros(size(Img));
        N = sqrt(size(Img,1)*size(Img,2));
        for z=1:size(Img,3)
            S=N./sum(sum(Img(:,:,z)));
            G(:,:,z) =double(Img(:,:,z)).*S - 1/N;
        end
    end
end
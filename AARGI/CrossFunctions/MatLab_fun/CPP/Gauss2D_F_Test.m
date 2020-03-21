%%
clear all;home;
[XI,YI] = meshgrid(1:10,1:10);
A=10.1;X=5.3;Y=6.2;Dx=1.5;Dy=3.5;Off=11.2;
Z=A*exp(-((XI-X)/Dx).^2/2-((YI-Y)/Dy).^2/2)+ Off;
Z2=Z.*(1+rand(10)*0.05);
surf(Z2);
%%
%%A = [9 4 1 5 1 10];
A = [5,5,5,5,5,5];
%[10.3395,5.2968,1.5054,6.2093,3.4664,11.5053;]
F = [1 1 1 1 1 1]; E = [0 0 0 0 0 0]; 
[Chi n] = GAUSS2D_vF2( [XI,YI] , Z2, A, F, E)
%%[10.4021,5.2936,1.4846,6.1797,3.4522,11.5091;]
%%
tifA=rand(256);
tifA_S=myShift1_vF(tifA, [-1 -1]);
tif1=myFFTNormal_h(tifA(11:138,11:138));
tif2=myFFTNormal_h(tifA_S(11:138,11:138));
tif1(256,256)=0;tif2(256,256)=0;
%%
tifShift = abs(fftshift(ifft2(fft2(tif1).*conj(fft2(tif2)))));
%tifShift = mycrossCorr2D(tif1,tif2);
[r c] = find(tifShift == max(tifShift(:)));
nHalf = length(tifShift)/2 + 1;
ss=-6:6;
tifShiftP=tifShift(ss+r(1),ss+c(1));
[XI,YI] = meshgrid(ss+7,ss+7);
A = [max(tifShiftP(:))-min(tifShiftP(:)),7,1,7,1,min(tifShiftP(:))];
F = [1 1 1 1 1 1]; E = [0 0 0 0 0 0];
[Chi n] = Gauss2D_vF2( [XI,YI] , tifShiftP, A, F, E);
fprintf('[%f %f]\n',(A(4)-7+(r-nHalf)), (A(2)-7+(c-nHalf)));
%% To compile
mex -O -output Gauss2D_vF2 Gauss2D_F.cpp
%% To compile with fft
%mex -output FIT_EXP1_vF2 FIT_EXP1_vF.cpp myFFT_Base.obj
%% Get data
% load('FitTest_dd');
load('ExampleTrace.mat'); dd = etrace;
%% Fit
Y=dd(:);
% Finding Fit Starting Values                                                
%    [ AMP,  T1, T2     offset,    X0, slope] 
X=[1:numel(Y)]'; 
[ii, jj]=max(Y); 
kk=min(Y); 
ll=find(Y>mean(Y));ll=numel(ll); 
A  = [ii-kk, ll, ll/10, kk, jj, 0];
% Fitting
F=[1 1 1 1 1 0]; E=F*0; [chi, N]=FIT_K_2EXP_SHIFT(X, Y, A, F ,E);
% plot result
figure;plot(X,dd);hold all; 
AA=A; X=[1:0.1:X(end)]; 
Yf = myEXP_Shift([AA(1) AA(2) 0 AA(5)],X)-myEXP_Shift([AA(1) AA(3) 0 AA(5)],X) + AA(4);
plot(X,Yf)
%%
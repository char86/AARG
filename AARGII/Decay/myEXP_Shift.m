function Out = myEXP_Shift(TT, X)
% Using:
% Out=myEXP_Shift(A, X);
% A = [Amp Tau Offset X0];
% Out = Amp exp (-(X-X0)./Tau) + Offset  or 0 if X<X0
Out = X.*0;
ii=find (X>=TT(4));
Out(ii)=TT(1)*exp(-(X(ii)-TT(4))/TT(2))+TT(3);
end

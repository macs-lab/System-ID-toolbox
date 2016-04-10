function [N_d,D_d]=cont2disc(N_c,D_c,Ts)
%function [N_d,D_d]=cont2disc(N_c,D_c,Ts)
%
%the function realize the conversion of a continues transfer function (in s):
%          bms^m + ... + b2s^2 + b1s + b0
%   G(s)= --------------------------------   ; m<=n
%         ans^n + ... + a2s^2 + a1s + a0
%to its 
%discret form (in z):
%          g0 + g1z^-1 + g2z^-2 + ... + gks^-k
%   G(z)= ------------------------------------   
%         h0 + h1z^-1 + h2z^-2 + ... + hls^-l
%using zero order hold.
%inputs:
%N_c = [bm ... b2 b1 b0] ... vector of continues numerator coefficients
%D_c = [an ... a2 a1 a0] ... vector of continues denominator coefficients
%Ts ... desired sampling time for discret transfer function
%outputs:
%N_d = [g0 g1 g2 ... gk] ... vector of discret numerator coefficients
%D_d = [h0 h1 h2 ... hl] ... vector of discret denominator coefficients
%
%written by: H. Prochazka, I.D. Landau
%7th june 2002

sys_c=tf(N_c,D_c);				%continue system
sys_d=c2d(sys_c,Ts,'zoh');		%conversion to a discret form (in z) with 0 order hold 
set(sys_d,'variable','z^-1');	%conversion to a discret form in z^-1
%-----
N_d=get(sys_d,'num');			%extraction and	
N_d=cat(1,N_d{:});				%conversion to obtain a vector of numerator coefficients 
D_d=get(sys_d,'den');			%extraction and	
D_d=cat(1,D_d{:});				%conversion to obtain a vector of denominator coefficients
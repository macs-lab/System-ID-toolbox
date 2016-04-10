function [re,im]=nyquist_ol(B,A,R,S,Ts);
%function [re,im]=nyquist_ol(B,A,R,S,Ts);
%plots a nyquist plot (frequency response in complex plain of open loop) 
%of a discret system consisted of the model B/A connected with a controller 
%R/S. The diagram is as follows:
%
% r     -----  u  -----   y
% --O--| 1/S |---| B/A |----
%  -|   -----     -----   | 
%   |   -----             |
%   ---|  R  |-------------
%       -----
%
%inputs:
%B=[b0 b1 ... bm] ... vector of model numerator coefficients (numerator: b0 + b1z^-1 + ...) 
%A=[a0 a1 ... an] ... vector of model denominator coefficients (denominator: a0 + a1z^-1 + ...)
%R=[r0 r1 ... rk] ... vector of controller numerator coefficients (numerator: r0 + r1z^-1 + ...)
%S=[s0 s1 ... sl] ... vector of controller denominator coefficients (denominator: s0 + s1z^-1 + ...)
%Ts ... sampling time of the system
%otputs:
%re ... vector of real parts of the frequency response
%im ... vector of imaginary parts of the frequency response
%
%written by:  H. Prochazka, I.D. Landau
%7th june 2002

num=conv(B,R);
den=conv(A,S);
sys=tf(num,den,Ts);
[re,im]=nyquist(sys);
nyquist(sys);
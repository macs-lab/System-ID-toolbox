function [num,den]=fd2pol(fnat,dpm,Ts);
%function [num,den]=fd2pol(fnat,dpm,Ts);
%computes numerator num and denominator den of 2nd order discret transfer function
%from natural frequency fnat and damping dmp of its continues equivalent.
%The continues 2nd order transfer function is of the form:
%                        2.pi.fnat
%   G2(s)= ---------------------------------------
%          s^2 + 2.dmp.2.pi.fnat.s + (2.pi.fnat)^2
%
%The discret transfer function is:
%            b1z^-1 + b2z^-2
%   Gd(z)= -------------------
%          1 + a1z^-1 + a2z^-1
%
%inputs:
%fnat ... natural frequency of continues 2nd order system
%dmp ... damping of continues 2nd order system
%Ts ... desired sampling time
%outputs:
%num=[0 b1 b2] ... vector of discret numerator coefficients
%den=[1 a1 a2] ... vector of discret denominator coefficients
%
%written by:  H. Prochazka, I.D. Landau
%7th june 2002


w=fnat*2*pi;
omega=w*sqrt(1-dpm^2);
alpha=exp(-dpm*w*Ts);
beta=cos(omega*Ts);
delta=sin(omega*Ts);
A1=-2*alpha*beta;
A2=alpha*alpha;
B1=1-alpha*(beta+(dpm*w*delta/omega));
B2=alpha^2+alpha*((dpm*w*delta/omega)-beta);
den=[1 A1 A2];
num=[0 B1 B2];
   
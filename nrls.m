function [B,A]=nrls(y,u,na,nb,d);
%function [B,A]=nrls(y,u,na,nb,d);
%linear digital model identification using non-recursive least square methode. The least square
%problem is resolved by QR-triangularization. The identified model is digital transfer 
%function of the following forme:
%                              z^-d (b_1z^-1 + b_2z^-2 + ... + b_nbz^-nb)
%                     G(z^-1)= ------------------------------------------
%                               1 + a_1z^-1 + a_2z^-2 + ... + a_naz^-na
%
%where na is the order of denominator, nb is the order of numerator and d is system's delay.
%inputs:
%y ...	vector of system outputs
%u ...	vector of system inputs
%na ...	order of model denominator
%nb ...	order of model numerator
%d ...	model's delay
%outputs:
%B ...	vector of numerator coefficients B=[0 0 ... 0 b_1 b_2 ... b_nb] (zeros is d+1)
%A ...	vector of denominator coefficients A=[1 a_1 a_2 ... a_na]
%
%written by:  H. Prochazka, I.D. Landau
%7th june 2002

sz=size(y);
if(sz(2)~=1), y=y';end;%column
sz=size(u);
if(sz(2)~=1), u=u';end;

nmax=max(na,nb+d);
ntheta=na+nb;
N=length(y)-nmax-1;

%y vector
Rnny=[];
for k=2:na+1,
   yp=y(k:k+N-1);
   Rnny=[yp Rnny];
end;
%u vector
Rnnu=[];
for k=d+2:nb+d+1,
   up=u(k:k+N-1);
   Rnnu=[up Rnnu];
end;
y0=y(1:N);   
Rnn=[Rnny Rnnu];
%QR-triangularization
[Q,Res]=qr(Rnn);
Y0=Q*y0;
Res1=Res(1:2*ntheta,:);
Y01=Y0(1:2*ntheta);
theta = Res\(Res'\(Rnn'*y0));
%refinement
ther=y0 - Rnn*theta;
er=Res\(Res'\(Rnn'*ther));
theta=theta + er;

A=[1 theta(1:na)'];
B=[zeros(1,d) 0 theta(na+1:na+nb)'];

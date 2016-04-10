function [F,E,EBst]=predisol(B,A,P);
%function [F,E,EBst]=predisol(B,A,P);
%solves the system of equations for predictor. Supposing a discret model of the form:
%   A(q^-1)y(t+d+1)=B*(q^-1)u(t)
%where 
%      B(q^-1)=q^-d B*(q^-1). 
%In general, this expression of filtred prediction is required:
%   P(q^-1)y_prdc(t+d+1)=fp[y(t),y(t-1),...,u(t),u(t-1),...]
%where P is a polynomial monique (le coeff. q^0 is 1) and stable. Then, the predictor has 
%following description:
%   P(q^-1)y(t+d+1)=F(q^-1)y(t)+E(q^-1)B*(q^-1)u(t)
%inputs:
%B=[ b0 b1 ... bm] ... vector of model numerator coefficients  (numerator: b0 + b1z^-1 + ...) 
%A=[a0 a1 ... an] ... vector of model denominator coefficients (denominator: a0 + a1z^-1 + ...) 
%P=[1 p2 p3 ... pk] ... polynomial filter of the predictor
%outputs:
%F=[f0 f1 f2 ... fl] ... vector of predictor coefficients (F(q^-1)=f0 + f1q^-1 + ...)
%E=[1 e1 e2 ... ew] ... vector of predictor coefficients (E(q^-1)=1 + e1q^-1 + ...)
%EBst=E(q^-1)B*(q^-1)
%
%written by:  H. Prochazka, I.D. Landau
%7th june 2002

PREC=1e-8;
%delay
d=1;
while abs(B(d))<PREC, d=d+1;end;
d=d-2;
%Bstar
Bst=B(d+2:length(B));
%removing zero at the end of vector
while abs(A(length(A)))<PREC, A=A(1:length(A)-1);end;
%M matrix construction
nA=length(A)-1;
nm=nA+d+1;
M=[];
for k=1:d+1,
   v=[zeros(k-1,1) ; A' ; zeros(d-k+1,1)];
   M=[M v];
end;
for k=1:nA,
   v=zeros(nm,1);
   v(d+1+k)=1;
   M=[M v];
end;
%adjusting the length of P
if length(P)<nm, P=[P PREC*ones(1,nm-length(P))]; end;
if length(P)>nm, P=P(1:nm);end;
P=P';
%inverse Mx=P => x=P/M
x=M\P;
E=x(1:d+1)';
F=x(d+2:length(x))';
EBst=conv(E,Bst);


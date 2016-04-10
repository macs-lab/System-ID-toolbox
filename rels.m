function [B,A,C]=rels(y,u,na,nb,nc,d,Fin,lam1,lam0)
% RELS  is used to identify a discrete time model of a plant operating in
%               open-loop based on the recursive extended least squares method.
%
%               [B,A,C]=rels(y,u,na,nb,nc,d,Fin,lam1,lam0)
%
%               y and u are the column vectors containing respectively the output and the
%               excitation signal.
%
%               na, nb are the order of the polynomials A,B and d is the pure time delay.
%               nc is the order of the numerator of the noise model.
%
%               Fin is the initial gain F0=Fin*(na+nb)*eye(na+nb) (Fin=1000 by default)
%
%               lam1 and lam0 make different adaptation algoritms as follows:
%
%               lam1=1;lam0=1           :decreasing gain (default algorithm)
%               0.95<lam1<1;lam0=1 :decreasing gain with fixed forgetting factor
%               0.95<lam1,lam0<1        :decreasing gain with variable forgetting factor
%
%               See also AFOLOE, FOLOE, OLOE, RLS and XOLOE.
%
%					 written by: A. Karimi, I.D. Landau
%					 7th june 2002



[nl,ncy]=size(y);
if ncy>2, error('This routine is only for SISO systems'),end
[nl,ncu]=size(u);
if ncu>2, error('This routine is only for SISO systems'),end
if (na<0 | nb<0 | nc<0 | d<0), error('The order of A,B,C and d should not be negative!!'),end


nd=min(length(u),length(y));    % number of data
nth=na+nb+nc;


if nd<nth, error('Number of data should be greater than the number of parameters!'),end



np=max([na+1,nb+d,nc+1]);


if nargin<6, error('This routin needs more parameters!'),end
if nargin<7, lam1=1;lam0=1;Fin=1000;end
if nargin<8, lam1=1;lam0=1;end
if nargin<9, lam0=1;end


if isempty(Fin), Fin=1000;end
if isempty(lam1),lam1=1;end
if isempty(lam0), lam0=1;end


if (lam1>1 | lam0>1), error('lam1 and lam0 should be less than 1');end
if (lam1<0.95 | lam0<0.95), disp ('warning :lam1 and lam0 are normally greater than 0.95');end


theta=zeros(nth,1);


phi=zeros(3*np,1);


F=Fin*nth*eye(nth);


i=[1:np-1 np+1:2*np-1 2*np+1:3*np-1];           %shift index for vector phi
j=[1:na np+d+1:np+d+nb 2*np+1:2*np+nc];                 %phi index for theta


for t=1:nd
        yhat=theta'*phi(j);
        e_apri=y(t)-yhat;
        e_apost=e_apri/(1+phi(j)'*F*phi(j));
        theta=theta+F*phi(j)*e_apost;


        F=1/lam1*(F-((F*phi(j)*phi(j)'*F)/(lam1+phi(j)'*F*phi(j))));
        lam1=lam0*lam1+1-lam0;


        yhat=theta'*phi(j);


        phi(i+1)=phi(i);phi(1)=-y(t);phi(np+1)=u(t);phi(2*np+1)=e_apost;


end
A=[1;theta(1:na)]';
B=[zeros(d+1,1);theta(na+1:na+nb)]';
C=[1;theta(na+nb+1:nth)]';

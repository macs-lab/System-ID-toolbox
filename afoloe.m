function varargout = afoloe(y,u,na,nb,d,Fin,lam1,lam0)
% AFOLOE        is used to identify a discrete time model of a plant operating in
%               open-loop based on the adaptive filtered output error method.
%
%               [B,A]=afoloe(y,u,na,nb,d,Fin,lam1,lam0)
%
%               y and u are the column vectors containing respectively the output and the
%               excitation signal.
%
%               na, nb are the order of the polynomials A,B and d is the pure time delay
%
%
%               Fin is the initial gain F0=Fin*(na+nb)*eye(na+nb) (Fin=1000 by default)
%
%               lam1 and lam0 make different adaptation algoritms as follows:
%
%               lam1=1;lam0=1           :decreasing gain (default algorithm)
%               0.95<lam1<1;lam0=1 :decreasing gain with fixed forgetting factor
%               0.95<lam1,lam0<1        :decreasing gain with variable forgetting factor
%
%               See also FOLOE, OLOE, XOLOE, RLS and RELS.
%
%					 written by: A. Karimi, I.D. Landau
%					 7th june 2002
% 
%                   Modified by Xu Chen
%                   mail.xchen@gmail.com
%                   2011-09-05

[nl,nc]=size(y);
if nc>2, error('This routine is only for SISO systems'),end
[nl,nc]=size(u);
if nc>2, error('This routine is only for SISO systems'),end
if (na<0 || nb<0 || d<0), error('The order of A,B and d should not be negative!!'),end

nd=min(length(u),length(y));    % number of data
nth=na+nb;

if nd<nth, error('Number of data should be greater than the number of parameters!'),end

np=max(na+1,nb+d);

if nargin<5, error('This routin needs more parameters!'),end
if nargin<6, lam1=1;lam0=1;Fin=1000;end
if nargin<7, lam1=1;lam0=1;end
if nargin<8, lam0=1;end

if isempty(Fin), Fin=1000;end
if isempty(lam1),lam1=1;end
if isempty(lam0), lam0=1;end

if (lam1>1 || lam0>1), error('lam1 and lam0 should be less than 1');end
if (lam1<0.95 || lam0<0.95), disp ('warning :lam1 and lam0 are normally greater than 0.95');end

theta=zeros(nth,1);
theta_vec = zeros(nth,nd);

phi=zeros(2*np,1);
phi_f=zeros(2*np,1);

F=Fin*nth*eye(nth);

i=[1:np-1 np+1:2*np-1];         %shift index for vector phi
j=[1:na np+d+1:np+d+nb];        %phi index for theta

for t=1:nd
    yhat=theta'*phi(j);
    e_apri=y(t)-yhat;
    e_apost=e_apri/(1+phi_f(j)'*F*phi_f(j));
    theta=theta+F*phi_f(j)*e_apost;
    
    theta_vec(:,t) = theta(:);
   
    F=1/lam1*(F-((F*phi_f(j)*phi_f(j)'*F)/(lam1+phi_f(j)'*F*phi_f(j))));
    lam1=lam0*lam1+1-lam0;
    
    yhat=theta'*phi(j);
    
    %Stability test for A
    A2=[1;theta(1:na)]';
    if max(abs(roots(A2))) < 1,A=A2;end
    
    % Filtering phi by 1/A
    phi(i+1)=phi(i);phi(1)=-yhat;phi(np+1)=u(t);
    yuf=[phi(1) phi(np+1)]-[A(2:na+1) zeros(1,np-na-1)]*[phi_f(1:np-1) phi_f(np+1:2*np-1)];
    phi_f(i+1)=phi_f(i);phi_f(1)=yuf(1);phi_f(np+1)=yuf(2);
end

A=[1;theta(1:na)]';
B=[zeros(d+1,1);theta(na+1:na+nb)]';

if nargout == 2
    varargout{1} = B;
    varargout{2} = A;
elseif nargout == 3
    varargout{1} = B;
    varargout{2} = A;
    varargout{3} = theta_vec;
else
    error('Incorrect number of outputs');
end
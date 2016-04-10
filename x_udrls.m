function varargout=x_udrls(y,u,na,nb,d,Fin,lam1,lam0)
% UDRLS   is used to identify a discrete time model of a plant operating in
% open-loop based on the recursive least squares method and the UD factorization
% for updating the adaption gain matrix (F).
%
%     [B,A]=udrls(y,u,na,nb,d,Fin,lam1,lam0)
%
%     y and u are the column vectors containing respectively the output and the
%     excitation signal.
%
%     na, nb are the order of the polynomials A,B and d is the pure time delay
%
%
%     Fin is the initial gain F0=Fin*(na+nb)*eye(na+nb) (Fin=1000 by default)
%
%     lam1 and lam0 make different adaptation algoritms as follows:
%
%     lam1=1;lam0=1           :decreasing gain (default algorithm)
%     0.95<lam1<1;lam0=1 :decreasing gain with fixed forgetting factor
%     0.95<lam1,lam0<1        :decreasing gain with variable forgetting factor
%
%     See also RLS, FOLOE, AFOLOE, XOLOE, OLOE and RELS.
%
%     written by: A. Karimi, I.D. Landau, F. Bouziani
%     7th june 2005
%     Updated by: Xu Chen mail.xchen@gmail.com 2012-03-27



[nl,nc]=size(y);
if nc>2, error('This routine is only for SISO systems'),end
[nl,nc]=size(u);
if nc>2, error('This routine is only for SISO systems'),end
if (na<0 | nb<0 | d<0), error('The order of A,B and d should not be negative!!'),end


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


if (lam1>1 | lam0>1), error('lam1 and lam0 should be less than 1');end
if (lam1<0.95 | lam0<0.95), disp ('warning :lam1 and lam0 are normally greater than 0.95');end

i=[1:np-1 np+1:2*np-1];         %shift index for vector phi
j=[1:na np+d+1:np+d+nb];        %phi index for theta

theta=zeros(nth,1);
phi=zeros(2*np,1);

U=eye(nth);
D=Fin*nth*eye(nth);
F=U*D*U';

theta_vec = zeros(nth,nd);

for t=1:nd
    yhat=theta'*phi(j);
    e_apri=y(t)-yhat;
    %     Beta(1)=1;
    Delta(1)=lam1;
    V=U'*phi(j);
    G=D*V;
    for J=1:nth
        %         Beta(J+1)=Beta(J)+ V(J)*G(J);
        Delta(J+1)=Delta(J)+ V(J)*G(J);
        D(J,J)=Delta(J)*D(J,J)/(Delta(J+1)*lam1);
        Gamma(J)=G(J);
        M(J)=-V(J)/Delta(J);
        if J>1
            for I=1:J-1
                Ut=U(I,J);
                U(I,J)=U(I,J)+Gamma(I)*M(J);
                Gamma(I)=Gamma(I)+Ut*Gamma(J);
            end
        end
    end
    for II=1:nth
        Gamma(II)=Gamma(II)/Delta(nth+1);
        %         Gamma(II)=Gamma(II)/Beta(nth+1);
    end
    e_apost=e_apri/Delta(nth+1);
    theta=theta+F*phi(j)*e_apost;  % or : theta=theta+Gamma'*e_apri;
    theta_vec(:,t) = theta(:);
    F=U*D*U';
    lam1=lam0*lam1+1-lam0;
    phi(i+1)=phi(i);phi(1)=-y(t);phi(np+1)=u(t);
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
end
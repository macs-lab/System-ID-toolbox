function[N0,nA0,nB0,d0]=estorderls(y,u,nmax,sel)
% Function to estimate the orders of a discrete time system.
% [N0,nA0,nB0,d0]=estorderls(y,u,nmax,sel)
% 
% Inputs:
% y = output data vector
% u = input data vector
% nmax = maximum system order (default 20)
% If sel is 1, the user will be asked to select the order N of the system 
% based on the order estimation curve.
% By default, sel=0. The value of N, minimizing the estimation criterion
% will be used.
% 
% Outputs:
% N0 = global order of the system
% nA0 = order of the denominator of the system
% nB0 = order of the numerator of the system
% d0 = identified delay
%
% Written by:  T.B. Airimitoaie, I.D. Landau and H. Prochazka
% 17th december 2009
% Last updated 15th february 2011

if nargin<2, error('"estorderlsg" needs a minimum of 2 parameters. Write "help estorderlsg" to the console for more informations.'); end
if nargin==2, nmax=20; sel=0; end
if nargin==3, sel=0; end

[V,S,VS]=estorderN(y,u,nmax);
figure
plot(0:nmax,VS,'r',0:nmax,V,'b')
title('Least Squares estimation of system order (N)')
if sel
    N0=input('Please input here the order of the system N = ');
else
    [vmin,imin]=min(VS);
    N0 = imin-1;
    disp(['N = ' int2str(N0)]);
end

[V,S,VS]=estdls(y,u,N0);
figure
plot(0:N0,VS,'r',0:N0,V,'b')
title('LS estimation of N-d')
[vmin2,imin2]=min(VS);
d0=N0-(imin2-1);
disp(['d = ' int2str(d0)]);

[V,S,VS]=estnBls(y,u,N0,d0);
figure
[vmin2,imin2]=min(VS);
nB0=imin2-1-d0;
plot(0:length(VS)-1,VS,'r',0:length(VS)-1,V,'b')
title('LS estimation of n_B+d')
disp(['nB = ' int2str(nB0)]);

[Va,Sa,VSa]=estnAls(y,u,N0,d0,nB0);
figure
plot(0:length(VSa)-1,VSa,'r',0:length(Va)-1,Va,'b')
title('LS estimation of n_A')
[vmin2,imin2]=min(VSa);
nA0 = imin2-1;
disp(['nA = ' int2str(nA0)]);

function [V,S,VS]=estorderN(y,u,nmax);
sz=size(y);
if(sz(2)~=1), y=y';end;%column
sz=size(u);
if(sz(2)~=1), u=u';end;

yz=[zeros(nmax-1,1);y];
uz=[zeros(nmax-1,1);u];

N=length(y)-nmax;
y0=y(nmax+1:N+nmax);   
V=sum(y0.^2)/N;
S=0;
VS=V;
for nn=1:nmax,
   %construction of R(nn)
   Rnny=[];
   Rnnu=[];
   for k=1:nn,
      yp=y(k:k+N-1);
      up=u(k:k+N-1);
      Rnny=[yp Rnny];
      Rnnu=[up Rnnu];
   end;
	y0=y(nn+1:N+nn);   
   Rnn=[Rnny Rnnu];
   theta_min=inv(Rnn'*Rnn)*Rnn'*y0;
   criter=((y0-Rnn*theta_min)'*(y0-Rnn*theta_min))/N;
   V=[V;criter];
   penal=0.4*nn*log10(N)/N;
   S=[S;penal];
end;
%normalization
V=V/V(1);
VS=V+S;
VS=VS/VS(1);

function [V,S,VS]=estdls(y,u,nmax)
sz=size(y);
if(sz(2)~=1), y=y';end;%column
sz=size(u);
if(sz(2)~=1), u=u';end;

yz=[zeros(nmax-1,1);y];
uz=[zeros(nmax-1,1);u];

N=length(y)-nmax;
V=[];
S=[];
VS=[];
for nn=0:nmax,
   %construction of R(nn)
   Rnny=[];
   Rnnu=[];
   for k=1:nmax,
      yp=y(k:k+N-1);
      Rnny=[yp Rnny];
   end;
   for k=1:nn,
      up=u(k:k+N-1);
      Rnnu=[up Rnnu];
   end;
	y0=y(nmax+1:N+nmax);   
   Rnn=[Rnny Rnnu];
   if size(Rnn)==[0 0], Rnn=1; end
   theta_min=inv(Rnn'*Rnn)*Rnn'*y0;
   criter=((y0-Rnn*theta_min)'*(y0-Rnn*theta_min))/N;
   V=[V;criter];
   penal=0.4*nn*log10(N)/N;
   S=[S;penal];
end;
%normalization
V=V/V(1);
VS=V+S;
VS=VS/VS(1);

function [V,S,VS]=estnBls(y,u,nmax,d0)
sz=size(y);
if(sz(2)~=1), y=y';end;%column
sz=size(u);
if(sz(2)~=1), u=u';end;

yz=[zeros(nmax-1,1);y];
uz=[zeros(nmax-1,1);u];

N=length(y)-nmax;
V=[];
S=[];
VS=[];
for nn=0:nmax,
   %construction of R(nn)
   Rnny=[];
   Rnnu=[];
   for k=1:nmax,
      yp=y(k:k+N-1);
      Rnny=[Rnny yp];
   end;
   for k=nmax:-1:nmax+1-nn%nmax-d0:-1:nmax-d0+1-nn,
      up=u(k:k+N-1);
      Rnnu=[Rnnu up];
   end;
	y0=y(nmax+1:N+nmax);   
   Rnn=[Rnny Rnnu];
   if size(Rnn)==[0 0], Rnn=1; end
   theta_min=inv(Rnn'*Rnn)*Rnn'*y0;
   criter=((y0-Rnn*theta_min)'*(y0-Rnn*theta_min))/N;
   V=[V;criter];
   penal=0.4*nn*log10(N)/N;
   S=[S;penal];
end;
%normalization
V=V/V(1);
VS=V+S;
VS=VS/VS(1);

function [V,S,VS]=estnAls(y,u,nmax,d0,nB0)
sz=size(y);
if(sz(2)~=1), y=y';end;%column
sz=size(u);
if(sz(2)~=1), u=u';end;

yz=[zeros(nmax-1,1);y];
uz=[zeros(nmax-1,1);u];

N=length(y)-nmax;
V=[];
S=[];
VS=[];
for nn=0:nmax,
   %construction of R(nn)
   Rnny=[];
   Rnnu=[];
   for k=nmax+1-nn:nmax,
      yp=y(k:k+N-1);
      Rnny=[yp Rnny];
   end;
   for k=nmax-d0:-1:nmax-d0-nB0+1
      up=u(k:k+N-1);
      Rnnu=[up Rnnu];
   end;
	y0=y(nmax+1:N+nmax);   
   Rnn=[Rnny Rnnu];
   if size(Rnn)==[0 0], Rnn=1; end
   aux = [y0 Rnn];
   aux(1:5,:);
   theta_min=inv(Rnn'*Rnn)*Rnn'*y0;
   criter=((y0-Rnn*theta_min)'*(y0-Rnn*theta_min))/N;
   V=[V;criter];
   penal=0.4*nn*log10(N)/N;
   S=[S;penal];
end;
%normalization
V=V/V(1);
VS=V+S;
% VS=VS/VS(1);
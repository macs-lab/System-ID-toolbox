function[N0,nA0,nB0,d0]=estorderiv(y,u,nmax,sel)
% Function to estimate the orders of a discrete time system.
% [N0,nA0,nB0,d0]=estorderiv(y,u,nmax,sel)
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

if nargin<2, error('"estorderivg" needs a minimum of 2 parameters. Write "help estorderivg" to the console for more informations.'); end
if nargin==2, nmax=20; sel=0; end
if nargin==3, sel=0; end

[V,S,VS]=estorderN(y,u,nmax);
[VS V];
figure
plot(0:nmax,VS,'r',0:nmax,V,'b')
title('Instrumental Variables estimation of system order (N)')
if sel
    N0=input('Please input here the order of the system N = ');
else
    [vmin,imin]=min(VS);
    N0 = imin-1;
    disp(['N = ' int2str(N0)]);
end

[V,S,VS]=estdiv(y,u,N0);
[VS V];
figure
plot(0:length(VS)-1,VS,'r',0:length(V)-1,V,'b')
title('IV estimation of N-d')
[vmin2,imin2]=min(VS);
d0 = N0-imin2+1;
disp(['d = ' int2str(d0)]);

[V,S,VS]=estnBiv(y,u,N0,d0);
[VS V];
figure
plot(0:length(VS)-1,VS,'r',0:length(V)-1,V,'b')
title('IV estimation of n_B+d')
[vmin2,imin2]=min(VS);
nB0 = imin2-1-d0;
disp(['nB = ' int2str(nB0)]);

[Va,Sa,VSa]=estnAiv(y,u,N0, d0, nB0);
[VSa Va];
figure
plot(0:length(VSa)-1,VSa,'r',0:length(Va)-1,Va,'b')
title('IV estimation of n_A')
[vmin2,imin2]=min(VSa);
nA0 = imin2-1;
disp(['nA = ' int2str(nA0)]);

function [VLs,S,VS]=estorderN(y,u,nmax)
if nargin==3,
elseif nargin==4, 
else return;
end;

sz=size(y);
if(sz(2)~=1), y=y';end;%column
sz=size(u);
if(sz(2)~=1), u=u';end;

VLs=[];
S=[];
VS=[];
L=2*nmax;
N=length(y)+1-L;
%construction of Z and y0
Z=[];
for k=1:2*nmax
    Z=[Z u(2*nmax-k+1:2*nmax-k+N-1)];
end
y0=y(L+1:length(y));

%QR-triangularization
[Q,R]=qr(Z);
Z1=R(1:L,:);
Zm=Z1'\Z';
Y0=Zm*y0;
for nn=0:nmax
    %construction of R(nn)
   Rnny=[];
   Rnnu=[];
   for k=1:nn
      yp=y(L-k+1:L-k+N-1);
      up=u(L-k+1:L-k+N-1);
      Rnny=[Rnny yp];
      Rnnu=[Rnnu up];
   end;
   if nn==0
       R1nn = 0;
       %theta computation
       theta_min=0;%inv(R1nn)*Y0;
   else
       R1nn=Zm*[Rnny Rnnu];
       %theta computation
       theta_min=(R1nn'*R1nn)\R1nn'*Y0;%inv(R1nn)*Y0;
   end
   criter=sum((Y0-R1nn*theta_min).^2)/N;
   VLs=[VLs;criter];
   penal=0.4*nn*log10(N)/N;
   S=[S;penal];
end
VLs=VLs/VLs(1);
VS=VLs+S;
VS=VS/VS(1);

function [VLs,S,VS]=estdiv(y,u,nmax)
if nargin==3,
elseif nargin==4, 
else return;
end;

sz=size(y);
if(sz(2)~=1), y=y';end;%column
sz=size(u);
if(sz(2)~=1), u=u';end;

VLs=[];
S=[];
VS=[];
L=2*nmax;
N=length(y)+1-L;   
%construction of R(nn) and Z
Z=[];
uz=u;
for k=1:2*nmax
    Z=[Z uz(2*nmax-k+1:2*nmax-k+N-1)];
end;
y0=y(L+1:length(y));

%QR-triangularization
[Q,R]=qr(Z);
Z1=R(1:L,:);
Zm=Z1'\Z';
Y0=Zm*y0;
for nn=0:nmax
   Rnny=[];
   Rnnu=[];
   for k=1:nmax
      yp=y(L-k+1:L-k+N-1);
      Rnny=[Rnny yp];
   end
   for k=1:nn
      up=u(L+(-nmax+k-1)+1:L+(-nmax+k-1)+N-1);
      Rnnu=[Rnnu up];
   end
   R1nn=Zm*[Rnny Rnnu];
   if size(R1nn)==[0 0], R1nn=1; end
	%theta computation
   theta_min=inv(R1nn'*R1nn)*R1nn'*Y0;%inv(R1nn)*Y0;
   criter=sum((Y0-R1nn*theta_min).^2)/N;
   VLs=[VLs;criter];
   penal=0.4*nn*log10(N)/N;
   S=[S;penal];
end
VLs=VLs/VLs(1);
VS=VLs+S;
VS=VS/VS(1);

function [VLs,S,VS]=estnBiv(y,u,nmax,d0)
if nargin~=4
    return;
end

sz=size(y);
if(sz(2)~=1), y=y';end;%column
sz=size(u);
if(sz(2)~=1), u=u';end;

VLs=[];
S=[];
VS=[];
L=2*nmax;
N=length(y)+1-L;   
%construction of R(nn) and Z
Z=[];
uz=u;
for k=1:2*nmax,
    Z=[Z uz(2*nmax-k+1:2*nmax-k+N-1)];
end;
y0=y(L+1:length(y));

%QR-triangularization
[Q,R]=qr(Z);
Z1=R(1:L,:);
Zm=Z1'\Z';
Y0=Zm*y0;

for nn=0:nmax
   Rnny=[];
   Rnnu=[];
   for k=1:nmax
      yp=y(L-k+1:L-k+N-1);
      Rnny=[Rnny yp];
   end
   for k=1:nn
      up=u(L-k+1:L-k+N-1);
      Rnnu=[Rnnu up];
   end
   R1nn=Zm*[Rnny Rnnu];
   if size(R1nn)==[0 0], R1nn=1; end
	%theta computation
   theta_min=inv(R1nn'*R1nn)*R1nn'*Y0;%inv(R1nn)*Y0;
   criter=sum((Y0-R1nn*theta_min).^2)/N;
   VLs=[VLs;criter];
   penal=0.4*nn*log10(N)/N;
   S=[S;penal];
end
VLs=VLs/VLs(1);
VS=VLs+S;
VS=VS/VS(1);

function [VLs,S,VS]=estnAiv(y,u,nmax,d0,nB0)
if nargin~=5
    return;
end

sz=size(y);
if(sz(2)~=1), y=y';end;%column
sz=size(u);
if(sz(2)~=1), u=u';end;

VLs=[];
S=[];
VS=[];
L=2*nmax;
N=length(y)+1-L;   
%construction of R(nn) and Z
Z=[];
for k=1:2*nmax,
    Z=[Z u(2*nmax-k+1:2*nmax-k+N-1)];
end;
y0=y(L+1:length(y));

%QR-triangularization
[Q,R]=qr(Z);
Z1=R(1:L,:);
Zm=Z1'\Z';
Y0=Zm*y0;

Rnnu=[];
for k=d0+1:d0+nB0
  up=u(L-k+1:L-k+N-1);
  Rnnu=[Rnnu up];
end

for nn=0:nmax
   Rnny=[];
   for k=1:nn
      yp=y(L-k+1:L-k+N-1);
      Rnny=[Rnny yp];
   end
   R1nn=Zm*[Rnny Rnnu];
   if size(R1nn)==[0 0], R1nn=1; end
   %theta computation
   theta_min=inv(R1nn'*R1nn)*R1nn'*Y0;%inv(R1nn)*Y0;
   criter=sum((Y0-R1nn*theta_min).^2)/N;
   VLs=[VLs;criter];
   penal=0.4*nn*log10(N)/N;
   S=[S;penal];
end
VLs=VLs/VLs(1);
VS=VLs+S;
VS=VS/VS(1);
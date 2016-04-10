function [wlossf,ulossf,wrni,urni,wyhat,uyhat]=olvalid(B,A,C,y,u)

%	OLVALID	is used for validation of plant models identified in open loop.
%			
%	[wlossf,ulossf,wrni,urni,wyhat,uyhat]=olvalid(B,A,C,y,u)
%
%   The estimation of the auto-correlation of the 
%	1-step prediction error and the estimated output. 
%
%	The estimation of the cross-correlation function between the 
%	output error prediction and the corresponding prediction error. 
%
%   The upper bound for the cross-correlation and autocorrelation estimations is given
%   for a confidence interval of %97.
%
%	The plant output y and the estimated ones (1-step predictor and output error predictor)
%   are compared and the loss functions are computed.
%			
%	B and A are respectively the numerator (including the pure time delay) and denominator of the plant 
%	model to be validated.
%   
%   C is the numerator of the noise model. If no noise model is avalaible C must be set to 1. 
%
%	y and u are the column vectors of output and input signals in
%	open loop operation.
%
%						
%	wyhat  : Open loop 1-step estimated output
%	wlossf : Loss function (the sum of (y-wyhat)^2 divided by number of data)
%
%	uyhat  : Open loop output error prediction
%	ulossf : Loss function (the sum of (y-uyhat)^2 divided by number of data)
%
%   wrni   : Auto-correlations
%   urni   : Uncorrelations 
%
%					 written by: G. Zito, I.D. Landau
%	                 Last modified on November 12, 2004 


na=length(A);
nB=length(B);

d=-1;i=1;
while B(i)==0,d=d+1;i=i+1;end
if d==-1, error('B should contain at least one zero!');end
if d==nB-1, error('B should be a nonzero vector!');end
nb=nB-d-1;


if isempty(C)
    nC=0;
else 
    nC=length(C);
end


nd=min(length(u),length(y));	% number of data
nc=max(max(nC,max(na-1,nb)),4);

    th_lim=num2str(2.17/sqrt(nd));

    disp(['   '])

    disp('Whiteness Test')    
    
    [wyhat,err_pred]=modelsim(B,A,C,y,u,d,nd);
    [wlossf,wrni]=correl(err_pred,err_pred,nc,nd)
    disp(['Theoretical upper limit : ' th_lim])
    disp(['   '])
    
    if (max(wrni)>0.15)
        beep;
        disp('WARNING!');
        disp('The model does not pass the whiteness test');
        disp(['   '])
    end
    
%%%%%%%%%%%%%%%%%%%%%

figure;
x1=1:nc;
xx=0;yy=0;
for i=1:nc,xx=[xx x1(i) x1(i) x1(i)];yy=[yy 0 wrni(i) 0];end
xx=[xx nc+1];yy=[yy 0];
plot(xx,yy,0:nc+1,2.17/sqrt(nd)*ones(1,length(0:nc+1)),0:nc+1,0.15*ones(1,length(0:nc+1)),'-')
hold on
xlabel('i');
ylabel('WRN(i)');
title('Whiteness test for the OL identification')
legend('Correlation terms',strcat('Theoretical limit: ',th_lim),'Practical limit: 0.15',0)

%%%%%%%%%%%%%%%%%%%%%
    
    disp(['   '])   
    disp('Uncorrelation Test')
    disp(['   '])
    
    sys_hat=filt(B,A,1);
    uyhat=lsim(sys_hat,u);
    err_dec=uyhat-y;
    ulossf=err_dec'*err_dec/nd
    [normr,urni]=correl(err_dec,uyhat,nc,nd);
    urni
    disp(['Theoretical upper limit : ' th_lim])
    disp(['   '])
    r0=abs(err_dec'*uyhat/nd);
    
    % Warning if  max(URNI) > 0.15 

    if (max(urni)>0.15)
        beep;
        disp('WARNING!');
        disp('The model does not pass the uncorrelation test');
        disp(['   '])
    end
    

   
%%%%%%%%%%%%%%%%%%%%%

figure;
x1=1:nc;
xx=0;yy=0;
for i=1:nc,xx=[xx x1(i) x1(i) x1(i)];yy=[yy 0 urni(i) 0];end
xx=[xx nc+1];yy=[yy 0];
plot(xx,yy,0:nc+1,2.17/sqrt(nd)*ones(1,length(0:nc+1)),0:nc+1,0.15*ones(1,length(0:nc+1)),'-')
hold on
xlabel('i');
ylabel('URN(i)');
title('Uncorrelation test for the OL identification')
legend('Correlation terms',strcat('Theoretical limit: ',th_lim),'Practical limit: 0.15',0)

%%%%%%%%%%%%%%%%%%%%%



% Warning if ( max(WRNI) and max(URNI) ) > 0.15 

if (max(wrni)>0.15 & max(urni)>0.15)
    beep;
    disp('WARNING!');
    disp('The model has to be improved');
end


function [y_pred,err_pred]=modelsim(B,A,C,y,u,d,nd)

na=length(A)-1;
nb=length(B)-1-d;


nc=length(C)-1;
nth=na+nb+nc;
np=max([na+1,nb+d,nc+1]);

i=[1:np-1 np+1:2*np-1 2*np+1:3*np-1];           %shift index for vector phi
j=[1:na np+d+1:np+d+nb 2*np+1:2*np+nc];                 %phi index for theta

theta=[A(2:end) B(2+d:end) C(2:end)]'; 
phi=zeros(3*np,1);

for t=1:nd
        y_pred(t)=theta'*phi(j);
        err_pred(t)=y(t)-y_pred(t);
        phi(i+1)=phi(i);phi(1)=-y(t);phi(np+1)=u(t);phi(2*np+1)=err_pred(t);
end

err_pred=err_pred';



function [normr,rni] = correl(x,y,nc,nd)

switch nargin
case 4
case 3
case 2
    nc=4;
case 1
    y=x;
    nc=4;
case 0
    error('This function needs at least one argument')
otherwise
    error('Too many arguments')
end

normr=sqrt((y'*y)*(x'*x))/nd;


for i=1:nc
    ri(i)=x(1+i:nd)'*y(1:nd-i)/(nd-i);
    rni(i)=abs(ri(i)/normr);
end


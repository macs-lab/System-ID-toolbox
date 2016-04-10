function [num,den,om_n,xi_n,om_d,xi_d]=filter22(Mt,ft,xid,Ts);
%function [num,den,om_n,xi_n,om_d,xi_d]=filter22(Mt,ft,xid,Ts);
%digital filter design function. Design bandstop digital filter
%with 2poles and 2zeros
%inputs:
%Mt ... desired attenuation (Mt<0)/amplification (Mt>0) in dB
%ft ... bandstop frequency in Hz
%xid ... damping of denominator poles for original analog filter. 
%        The resulting denominator damping of digital filter will
%        be in general close to this value. Recommended values 
%        are from 0.6 to 0.9
%Ts ... sampling time
%outputs:
%num ... numerator of digital filter
%den ... denominator of digital filter
%om_n ... natural frequency of numerator [rad/s]
%xi_n ... damping of numerator
%om_d ... natural frequency of denominator [rad/s]
%xi_d ... damping of denominator
%
%written by:  H. Prochazka, I.D. Landau
%7th june 2002

om_zt=2*pi*Ts*ft;

om_st=2*tan(om_zt/2)/Ts;

Mlin=10^(Mt/20);

xin=xid*Mlin;

bz0=4/Ts^2 + 4*xin*om_st/Ts + om_st^2;
bz1=2*om_st^2 - 8/Ts^2;
bz2=4/Ts^2 - 4*xin*om_st/Ts + om_st^2;

az0=4/Ts^2 + 4*xid*om_st/Ts + om_st^2;
az1=2*om_st^2 - 8/Ts^2;
az2=4/Ts^2 - 4*xid*om_st/Ts + om_st^2;

num=[bz0 bz1 bz2];
den=[az0 az1 az2];
%polynomials to frequency and damping
%zeros
zd1 = -num(2) + sqrt( num(2)^2-4*num(1)*num(3) );
zd1=zd1/( 2*num(1) );
Absd=abs(zd1);
Phsd=angle(zd1);
om_n=sqrt( (Phsd^2+log(Absd)^2)/Ts^2 );
xi_n=-log(Absd)/(om_n*Ts);
disp('ZEROS:');
disp(sprintf('nat.frequency[rad/s]\t nat.frequency[Hz]\t damping'));
disp(sprintf('%f\t\t %f\t\t %f\n',om_n,om_n/(2*pi),xi_n));
%poles
zd1 = -den(2) + sqrt( den(2)^2-4*den(1)*den(3) );
zd1=zd1/( 2*den(1) );
Absd=abs(zd1);
Phsd=angle(zd1);
om_d=sqrt( (Phsd^2+log(Absd)^2)/Ts^2 );
xi_d=-log(Absd)/(om_d*Ts);
disp('POLES:');
disp(sprintf('nat.frequency[rad/s]\t nat.frequency[Hz]\t damping'));
disp(sprintf('%f\t\t %f\t\t %f',om_d,om_d/(2*pi),xi_d));
%normalization
num=num/den(1);
den=den/den(1);

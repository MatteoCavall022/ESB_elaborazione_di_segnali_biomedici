function [teta] = delay(fft1r,fft1i,fft2r,fft2i,start)

% function [teta] = Delay(fft1r,fft1i,fft2r,fft2i,start)
% ff1 = ff1r + j ff1i
% ff2 = ff2r + j ff2i
%
% Computation of the delay between two waveforms using the spectrum alignment.
% Newton's method for the minimization of the least square criterium on FFT. 

 % initialization of the error plot

 cpter=1;

 % initialisation
 t=start;
 teta=10;
 n=length(fft1r);
 coef=2*pi/n;

 % loop of computation of the derives (de1 & de2)
 % of the alignment error while the improvment step
 % of teta is more than 0.05 
 while (abs(teta-t)>=5e-5)
  teta=t;
  de1=0;
  de2=0;
  b=coef*teta;
  for k=2:(n/2)
    cs=cos(2*(k-1)*pi*teta/n);	
    sn=sin(2*(k-1)*pi*teta/n);	
    f1tr=(fft1r(k)*cs-fft1i(k)*sn);
    f1ti=(fft1r(k)*sn+fft1i(k)*cs);
    de1=de1+(2*(k-1)*pi/n)*(f1tr*(-fft2i(k))+f1ti*fft2r(k));
    de2=de2+(2*(k-1)*pi/n)^2*(f1tr*fft2r(k)+f1ti*fft2i(k));
  end
  de1=de1*4/n;
  de2=de2*4/n;
  % Newton's criteria
  if (de2>0) 
    u=-de1/de2;
    if (abs(u)>0.5)
      u=-0.5*abs(de1)/de1;
    end
  else
    u=-0.5*abs(de1)/de1;
  end
  err1(cpter,1)=de1;
  err2(cpter,1)=de2;
  u_v(cpter,1)=u;
  cpter=cpter+1; 
  t=teta+u;	  % result
 end
% teta=t;
	

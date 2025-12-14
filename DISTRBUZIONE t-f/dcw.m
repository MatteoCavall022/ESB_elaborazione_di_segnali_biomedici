function w=dcw(x,fc,ntlag,sig,T);
% DCW - Discrete Choi-Williams time-frequency transform
% Usage:   w=dcw(x,fc,ntlag,sig,T);
% 
% Inputs:  x      input signal
%          ntlag  number of time lag for the autocorrelation sequence estimate
%          fc     sampling rate in Hz
%          sig    value of the selectivity parameter (>0 and <100)
%          T      time support of the input signal in seconds
% Output:  w      matrix containing the Choi-Williams transform

%	author(s):	G. Gagliati
%			P. Bonato, 4-22-1996, revised

aa=daf(x,ntlag,fc,T);
aa=1/((1/fc)*(1/T))*aa;
[rowaa colaa]=size(aa);
kcw=kercw(sig,rowaa,colaa,fc);
kcw=fftshift(kcw);
ww=kcw.*aa;
clear aa
clear kcw
ww=ww';
w=ifft(ww);
clear ww
fprintf('First Fourier transform:done!\n');
w=w';
www=1/(ntlag*2/fc)*(2/fc)*fft(w);
clear w
w=www;
clear www
fprintf('Second Fourier transform:done!\n');

[roww colw]=size(w);
w=w(1:roww/2+1,:);
fprintf('Choi-Williams estimate:done!\n');
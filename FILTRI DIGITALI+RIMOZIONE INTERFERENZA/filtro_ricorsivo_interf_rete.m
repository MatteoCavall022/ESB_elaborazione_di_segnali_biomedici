% filtro ricorsivo reiezione 50 Hz

function [b,a]=filtro_ricorsivo_interf_rete(z,B,fc,fs,ft,x,xor)

% input parameters
%	z  (0.01) attenuazione minima alla frequenza fc
%	B  (2-8) larghezza di banda corrispondente alla attenuazione 0.707
% 	fc (50) frequenza di centro banda
%fs, frequenza campionamento
%ft,frequenza di taglio 
%xor segnale pulito, x segnale sporco
% output parameters
%	cMA  filter coefficients (MA part)
%	cAR  filter coefficients (AR part)

k=fs/ft;
if k==floor(k)
    fs2=fs;
    xx=x;
    xxor=xor;
else 
    fs2=(floor(k)*ft);
    xx=resample(x,fs2,fs); %ricampioni a 500
    xxor=resample(xor,fs2,fs);
    subplot(2,1,1)
    plot(x(1:1024),'r')
    axis([0 1024 -300 500])
    title('Original ECG signal - fs= %d Hz', fs)
    subplot(2,1,2)
    plot(xx(1:1000),'g')
    axis([0 1000 -300 500])
    title('Resampled ECG signal - fs= %d Hz',fs2)

end
k=fs2/ft;
T=1/fs2; %intervallo di campionamento
b=pi*B*T;
a=b*z;
c1=-2*(1-a)*cos(2*pi*fc*T);
c2=(1-a)^2;
c3=2*(1-b)*cos(2*pi*fc*T);
c4=-(1-b)^2;
b=[1 c1 c2];
a=[1 -c3 -c4];
 figure(2)
 freqz(b,a,fs,fs2);
 figure(3)
 y=filtfilt(b,a,xx);
 subplot(3,1,1);
 plot(xxor(1:1500),'r');
 axis([0 1500 -500 500]);
 title('Original ECG signal (uncorrupted)')
 subplot(3,1,2);
 plot(x(1:1500),'y');
 title('Corrupted ECG signal (50 Hz interference)')
 subplot(3,1,3);
 plot(y(1:1500),'g');
 title('Filtered ECG signal (Recursive filter)')

end
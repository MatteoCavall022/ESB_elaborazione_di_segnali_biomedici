function [b,fs2,k] = notch_antisimmetrico (fs,ft,x,xor)
%fs, freq campionamento
%ft, freq taglio
%x segnale sporco, xor segnale pulito
k=fs/ft;
if k==floor(k)
    fs2=fs;
    x=xx;
    xor=xxor;
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
a=1;
b=[1, zeros(1,k-1), -1]; %fs/ft+1=11, gli zeri sono nove perchè devono essere k-1, k=500/50
figure (2)
freqz(b,a,fs,fs2)
 figure(3)
 y=filter(b,a,xx);
 subplot(3,1,1);
 plot(xxor(1:1500),'r'); %segnale pulito 
 axis([0 1500 -500 500]);
 title('Original ECG signal (uncorrupted)')
 subplot(3,1,2);
 plot(x(1:1500),'y'); %segnale sporco 
 title('Corrupted ECG signal (50 Hz interference)')
 subplot(3,1,3);
 plot(y(1:1500),'g'); %segnale filtrato con notch
 title('Filtered ECG signal (simple notch)')
end
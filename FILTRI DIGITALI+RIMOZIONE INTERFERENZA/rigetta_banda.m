function [b,a,fs2,y] = rigetta_banda(ord_filtro,fs, type_filter,ft,x,xor)
% ord_filtro, ordine filtro o grado o num coefficienti-1, da testo
%fs, f campionamento
%type_filter, FIR o IIR_doppia_passata
%ft, f taglio, per interferenza di rete=50
%se non sai cosa inserire guarda pdf LAB02 ultime righe
inizio_banda_pass=input('inserire inizio banda passante: ');
inizio_banda_attenuata=input('inserire inizio banda attenuata: ');
fine_banda_attenuata=input('inserire fine banda attenuata: ');
fine_banda_pass=input('inserire fine banda passante: ');
k=fs/ft;
if k==floor(k)
    fs2=fs;
    x=xx;
    xor=xxor;
else 
    fs2=(floor(k)*ft);
    xx=resample(x,fs2,fs); %ricampioni a 500
    xxor=resample(xor,fs2,fs);
    figure(1)
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
fNy=fs2/2;
if strcmp(type_filter, 'IIR_doppia_passata')
     f=[0, inizio_banda_pass, inizio_banda_attenuata, fine_banda_attenuata, fine_banda_pass, fNy]/fNy; 
     m=[1 1 0 0 1 1];
     [b,a]=yulewalk(ord_filtro,f,m);
     figure(2)
     freqz(b,a,fs,fs2);
     y=filtfilt(b,a,xx);
     figure(3)
     subplot(3,1,1);
     plot(xxor(1:1500),'r');
     axis([0 1500 -500 500]);
     title('Original ECG signal (uncorrupted)')
     subplot(3,1,2);
     plot(x(1:1500),'y');
     title('Corrupted ECG signal (50 Hz interference)')
     subplot(3,1,3);
     plot(y(1:1500),'g');
     title('Filtered ECG signal (Anticausal IIR band reject)') 

    %grafici sovrapposti
    figure(4)
    subplot(1,1,1);
    plot(xxor(1:500),'r')
    hold on
    plot(y(1:500),'y')
    title('Comparison between original and filtered signals IIR') 
    hold off
elseif strcmp(type_filter, 'FIR')
    f=[0 inizio_banda_pass inizio_banda_attenuata fine_banda_attenuata fine_banda_pass fNy]/fNy;
    m=[1 1 0 0 1 1];
    [b,a]=firpm(ord_filtro , f, m); 
    figure(2)
    freqz(b,a,fs,fs2);
    y=filter(b,a,xx);
    figure(3)
    subplot(3,1,1);
    plot(xxor(1:1500),'r');
    axis([0 1500 -500 500]);
    title('Original ECG signal (uncorrupted)')
    subplot(3,1,2);
    plot(x(1:1500),'y');
    title('Corrupted ECG signal (50 Hz interference)')
    subplot(3,1,3);
    plot(y(1:1500),'g');
    title('Filtered ECG signal (FIR band reject)')
    figure(4)
    %sovrapposizione grafici
    subplot(1,1,1);
    plot(xxor(1:500),'r')
    hold on
    plot(y(1:500),'y')
    title('Comparison between original and filtered signals (FIR band reject)') 
    hold off

end
end
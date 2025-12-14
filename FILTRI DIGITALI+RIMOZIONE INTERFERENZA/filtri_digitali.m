function [a,b,ytm] = filtri_digitali(xor,xtm,fs, typefilter)
%tramite load dei dati devi arrivare ad avere xor segnale pulito e xtm
%segnale corrotto da rumore- un es è in lab02cm_emg
%typefilter è tipo di filtro: IIR_doppia_passata, FIR,
%filtraggio_blocchi_FIR, filtraggio_blocchi_IIR
%in base al filtro che ti serve al posto di typefilter scrivi uno dei 3
%sopra
 zoom=1000; %numero indicativo per visualizzare un pezzo di segnale (qui avevo 5120 campioni e ho fatto 1000)
 disp ('usa il grafico per trovare il range in cui hai solo rumore e non segnale per darlo in input per PSD rumore')
 subplot(2,1,1)
 plot(xor(1:zoom),'r'); 
 axis([0 zoom -300 500]);
 title('Segnale ECG originale')
 subplot(2,1,2)
 plot(xtm(1:zoom),'g');
 axis([0 zoom -300 500]);
 title('Segnale ECG corrotto da tremore')
%STIMA DELLA PSD DEL RUMORE (scegli come f di taglio da dare tra poco in
%input dove hai il picco)

inizio_rumore=input('inserire punto inizio rumore > di 0: ');
fine_rumore=input('inserire punto fine rumore: ');
%hann (lunghezza finestra), overlap, NFFT
[Pxx,f] = pwelch((xtm(inizio_rumore:fine_rumore)-mean(xtm(inizio_rumore:fine_rumore))),hann(128),64,128,fs);
subplot(1,1,1) 
 plot(f,Pxx)
 axis([0 256 0 max(abs(Pxx))])
 title('Stima della PSD del rumore') 

 fNy=fs/2;
 ft=input('inserisci ft da Psd del rumore: '); %se hai ECG che arriva fino a 40, tagli a 50

if strcmp(typefilter, 'IIR_doppia_passata')
    funzione=input('inserisci il tipo di funz (yulewalk,cheby1,butter: ','s');
    if strcmp (funzione, 'yulewalk')
        f=[0,ft/fNy,(ft+10)/fNy,1];       %se passa-alto f=, se passabanda f=[0,ft,ft+10,ft+20,ft,fNy]/fNy               
        m=[1,1,0,0];                      %se filtro passa-alto m=[0,0,1,1] se passabanda m=[1,1,0,0,1,1]
        [b,a]=yulewalk(11,f,m);                  
        freqz(b,a,fs,fs)  % plot di modulo e fase del filtro progettato
        ytm=filtfilt(b,a,xtm);
    
    elseif strcmp (funzione, 'cheby1')
        ord_filtro=input('inserisci ordine filtro o (num_coefficienti-1): '); % se non dice niente, sui 5 circa
        tipo_filtraggio=input('bandpass,low,high: ','s');
        [b,a] = cheby1(ord_filtro, 0.5, [1, ft]/fNy,tipo_filtraggio);            
        freqz(b,a,fs,fs)
        ytm=filtfilt(b,a,xtm);
   
    elseif strcmp(funzione, 'butter')
        ord_filtro=input('inserisci ordine filtro o (num_coefficienti-1): ');
        tipo_filtraggio=input('low,high: ','s');
        [b, a] = butter(ord_filtro, ft/fNy, tipo_filtraggio);
        freqz(b,a,fs,fs)
        ytm=filtfilt(b,a,xtm);
end
  
 subplot(2,1,1) 
 plot(xor(1:zoom),'r') %plotto xor per capire se il mio filtraggio che ha prodotto ytm è vicino al segnale teorico pulito
 title('Segnale ECG originale')
 subplot(2,1,2)
 plot(ytm(1:zoom),'g')
 title ('Segnale filtrato (IIR  anticausale)')

elseif strcmp(typefilter, 'FIR')
    funzione=input('inserisci il tipo di funz (fir1,firpm: ','s');
    if strcmp(funzione, 'fir1')
        ord_filtro=input('inserisci ordine filtro o (num_coefficienti-1): ');
        tipo_filtraggio=input('bandpass,low,high: ','s');
        b=fir1(ord_filtro,[1,ft]/fNy, tipo_filtraggio);                
        a=1;
        freqz(b,a,fs,fs)
        ytm=filter(b,a,xtm); 
 
    elseif strcmp(funzione,'firpm')
        ord_filtro=input('inserisci ordine filtro o (num_coefficienti-1): ');
        F=[0 (ft/fNy) ((ft/fNy)+0.1) 1];
        A=[1 1 0 0]; %se passa alto [0 0 1]
        b = firpm(ord_filtro,F,A);  
        a=1;
        freqz(b,a,fs,fs)
        ytm=filter(b,a,xtm);
    end

                 
 subplot(2,1,1) 
 plot(xor(1:zoom),'r')
 axis([0 zoom -300 500])  
 title('Segnale ECG originale')
 subplot(2,1,2)
 plot(ytm(1:zoom),'g')
 axis([0 zoom -300 500]) 
 title ('Segnale filtrato')
 %ritardo introdotto dal filtro FIR è ord_filtro/2

elseif strcmp(typefilter, 'filtraggio_blocchi_FIR')
   ord_filtro=input('inserisci ordine filtro o (num_coefficienti-1): ');
   length_blocco=input('inserisci numero di campioni per blocco contiguo: ');
   ytm=(zeros(1,length(xtm)));     
   zi=(zeros(1,ord_filtro));
   blocco=length(xtm)/length_blocco;
   tipo_filtraggio=input('bandpass,low,high: ','s');
   b=fir1(ord_filtro,[1,ft]/fNy, tipo_filtraggio);
   a=1;
 for i=1:(blocco) %tutto segnale è 5120, dice di filtrare blocchi contigui di 512 campioni, quindi 5120/512=10
   [yy,zi]=filter(b,a,xtm(((i-1)*length_blocco+1):(i*length_blocco)), zi);                
   ytm(((i-1)*length_blocco+1):(i*length_blocco))=yy; 
 end
 subplot(2,1,1) 
 plot(xor(1:(3*length_blocco)),'r')
 title('Segnale ECG originale')
 subplot(2,1,2) 
 plot(ytm(1:1400),'g')
 title('Segnale filtrato')


elseif strcmp(typefilter, 'filtraggio_blocchi_IIR')
    ord_filtro = input('inserisci ordine filtro o (num_coefficienti-1): ');
    length_blocco = input('inserisci numero di campioni per blocco contiguo: ');
    ytm = zeros(1, length(xtm));     
    blocco = length(xtm) / length_blocco;
    tipo_filtraggio = input('bandpass,low,high: ','s');

    [b, a] = cheby1(ord_filtro, 0.5, [1, ft]/fNy, tipo_filtraggio); % togli 1 se metti low o high, [1, ft]/fNy è Wp
    nzi = max(length(a), length(b)) - 1;
    zi = zeros(1, nzi);

    for i = 1:blocco
        idx_start = (i-1)*length_blocco + 1;
        idx_end = i*length_blocco;
        seg = xtm(idx_start:idx_end);

        % Prima passata (diretta)
        [yy, zi1] = filter(b, a, seg, zi);

        % Seconda passata (inversa)
        yy_rev = flip(yy); % inverti il vettore
        [yy2, zi2] = filter(b, a, yy_rev, zi);
        yy_final = flip(yy2); % reinverti il risultato

        ytm(idx_start:idx_end) = yy_final;
        % opzionale: aggiorna zi se vuoi continuità tra blocchi (ma per filtfilt di solito non serve)
    end

    subplot(2,1,1) 
    plot(xor(1:(3*length_blocco)),'r')
    title('Segnale ECG originale')
    subplot(2,1,2) 
    plot(ytm(1:1400),'g')
    title('Segnale filtrato (doppia passata)')
end
end



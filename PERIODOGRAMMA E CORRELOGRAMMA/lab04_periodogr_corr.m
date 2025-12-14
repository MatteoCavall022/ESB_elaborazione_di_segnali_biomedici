clear all
close all
clc
%% Parametri del segnale
fs = 512;             % Frequenza di campionamento (Hz)
T = 15;               % Durata acquisizione (s)
nCanali = 25;         % Numero di canali
n_campioni = fs * T;           % Numero di campioni per ciascun canale
t=1;
N=fs*t;

%% Lettura del file binario multiplexato
fid = fopen('eeg2501.bin','rb');
if fid == -1
    error('File eeg2501.bin non trovato!');
end

% Leggi tutti i dati come float (single precision)
data = fread(fid, [nCanali, N], 'float32'); % [canali x campioni]
fclose(fid);

% Ogni riga di "data" è un canale, ogni colonna un campione

%% Estrazione dei canali 7 e 9
canale7 = data(7, :); % Vettore riga, campioni del canale 7
canale9 = data(9, :); % Vettore riga, campioni del canale 9

%% Estrazione del decimo secondo (da 9 a 10 s)
t_inizio = 9 * fs + 1;      % Indice di inizio (MATLAB parte da 1)
t_fine   = 10 * fs;         % Indice di fine

seg7 = canale7(t_inizio:t_fine);
seg9 = canale9(t_inizio:t_fine);
time = (0:length(seg7)-1) / fs + 9; % Asse temporale da 9 a 10 s

%% Visualizzazione
figure;
subplot(2,1,1)
plot(time, seg7)
xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title('EEG Canale 7 (decimo secondo)')

subplot(2,1,2)
plot(time, seg9)
xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title('EEG Canale 9 (decimo secondo)')

%% calcolo mediante periodogramma semplice
% se non ti dice niente NFFT=N altrimenti NFFT=fc/ris_app
NFFT=N;
[Pxx_7,f_7]=pwelch(canale7-mean(canale7),rectwin(N),0,NFFT,fs);
[Pxx_9,f_9]=pwelch(canale9-mean(canale9),rectwin(N),0,NFFT,fs);
subplot(2,1,1)
plot(f_7,Pxx_7 ./max(Pxx_7))
title ('Periodogramma semplice canale 7')
subplot(2,1,2)
plot(f_9, Pxx_9 ./max(Pxx_9))
title('periodogramma semplice canale 9')


%% calcolo della ripartizione della potenza media del segnale nelle bande δ (0, 5 Hz → 3, 5 Hz), θ (3, 5 Hz → 7 Hz), α (7 Hz → 14 Hz), β1 (14 Hz → 21 Hz) e β2 (21 Hz → 30 Hz).
%segnale 7
bande = [0.5, 3.5; 3.5, 7; 7, 14; 14, 21; 21, 30];
indici_bande = zeros(size(bande));
for i = 1:size(bande,1)
    [~, indici_bande(i,1)] = min(abs(f_7 - bande(i,1))); % Trova indice più vicino
    [~, indici_bande(i,2)] = min(abs(f_7 - bande(i,2)));
end
indici_bande = min(max(indici_bande, 1), length(Pxx_7));
disp(indici_bande)

 Ptot_7=sum(Pxx_7);
    P_7=sum(Pxx_7(indici_bande(1,1):indici_bande(1,2))); 
    Pd_7=P_7/Ptot_7;
    P_7=sum(Pxx_7(indici_bande(2,1):indici_bande(2,2)));
    Pt_7=P_7/Ptot_7;
    P_7=sum(Pxx_7(indici_bande(3,1):indici_bande(3,2)));
    Pa_7=P_7/Ptot_7;
    P_7=sum(Pxx_7(indici_bande(4,1):indici_bande(4,2)));
    Pb1_7=P_7/Ptot_7;
    P_7=sum(Pxx_7(indici_bande(5,1):indici_bande(5,2)));
    Pb2_7=P_7/Ptot_7; 
    fprintf('\n Pd_7 = %.4f Pt_7 = %f Pa_7 = %f Pb1_7 = %f Pb2_7 = %f \n',Pd_7,Pt_7,Pa_7,Pb1_7,Pb2_7);

    %segnale 9
    bande = [0.5, 3.5; 3.5, 7; 7, 14; 14, 21; 21, 30];
indici_bande = zeros(size(bande));
for i = 1:size(bande,1)
    [~, indici_bande(i,1)] = min(abs(f_9 - bande(i,1))); % Trova indice più vicino
    [~, indici_bande(i,2)] = min(abs(f_9 - bande(i,2)));
end
indici_bande = min(max(indici_bande, 1), length(Pxx_9));
disp(indici_bande)
    Ptot_9=sum(Pxx_9);
    P_9=sum(Pxx_9(indici_bande(1,1):indici_bande(1,2))); 
    Pd_9=P_9/Ptot_9;
    P_9=sum(Pxx_9(indici_bande(2,1):indici_bande(2,2)));
    Pt_9=P_9/Ptot_9;
    P_9=sum(Pxx_9(indici_bande(3,1):indici_bande(3,2)));
    Pa_9=P_9/Ptot_9;
    P_9=sum(Pxx_9(indici_bande(4,1):indici_bande(4,2)));
    Pb1_9=P_9/Ptot_9;
    P_9=sum(Pxx_9(indici_bande(5,1):indici_bande(5,2)));
    Pb2_9=P_9/Ptot_9; 
    fprintf('\n Pd_9 = %.4f Pt_9 = %f Pa_9 = %f Pb1_9 = %f Pb2_9 = %f \n',Pd_9,Pt_9,Pa_9,Pb1_9,Pb2_9);

    %% periodogramma con metodo di welch
    %metodo di Welch con brani della lunghezza di 1024 campioni sovrapposti di 512 campioni ed utilizzando
    % una finestra di Hamming (di lunghezza pari al brano). %qui
    % specificato lunghezza del brano in hammming altrimenti circa N/4
    N_welch=1024; % se non ti dice niente = N iniziale
    overlap=512; %solitamente del 50%
    NFFT=N_welch; %cioè uguale a N o N_welch in base ai casi sopra se non dice niente, altrimenti NFFT=fc/ris_app
    [Pxx_7,f_7]=pwelch(canale7-mean(canale7),hamming(N_welch),overlap,NFFT,fs);
    [Pxx_9,f_9]=pwelch(canale9-mean(canale9),hamming(N_welch),overlap,NFFT,fs);
    figure(2)
    subplot(2,1,1)
    plot(f_7,Pxx_7 ./max(Pxx_7))
    title ('Metodo di welch Canale 7') 
    subplot(2,1,2)
    plot(f_9,Pxx_9 ./max(Pxx_9))
    title ('Metodo di welch Canale 9')

    %% calcolo della ripartizione della potenza media del segnale nelle bande δ (0, 5 Hz → 3, 5 Hz), θ (3, 5 Hz → 7 Hz), α (7 Hz → 14 Hz), β1 (14 Hz → 21 Hz) e β2 (21 Hz → 30 Hz).
%segnale 7
bande = [0.5, 3.5; 3.5, 7; 7, 14; 14, 21; 21, 30];
indici_bande = zeros(size(bande));
for i = 1:size(bande,1)
    [~, indici_bande(i,1)] = min(abs(f_7 - bande(i,1))); % Trova indice più vicino
    [~, indici_bande(i,2)] = min(abs(f_7 - bande(i,2)));
end
indici_bande = min(max(indici_bande, 1), length(Pxx_7));
disp(indici_bande)

 Ptot_7=sum(Pxx_7);
    P_7=sum(Pxx_7(indici_bande(1,1):indici_bande(1,2))); 
    Pd_7=P_7/Ptot_7;
    P_7=sum(Pxx_7(indici_bande(2,1):indici_bande(2,2)));
    Pt_7=P_7/Ptot_7;
    P_7=sum(Pxx_7(indici_bande(3,1):indici_bande(3,2)));
    Pa_7=P_7/Ptot_7;
    P_7=sum(Pxx_7(indici_bande(4,1):indici_bande(4,2)));
    Pb1_7=P_7/Ptot_7;
    P_7=sum(Pxx_7(indici_bande(5,1):indici_bande(5,2)));
    Pb2_7=P_7/Ptot_7; 
    fprintf('\n Pd_7 = %.4f Pt_7 = %f Pa_7 = %f Pb1_7 = %f Pb2_7 = %f \n',Pd_7,Pt_7,Pa_7,Pb1_7,Pb2_7);

    %segnale 9
    bande = [0.5, 3.5; 3.5, 7; 7, 14; 14, 21; 21, 30];
indici_bande = zeros(size(bande));
for i = 1:size(bande,1)
    [~, indici_bande(i,1)] = min(abs(f_9 - bande(i,1))); % Trova indice più vicino
    [~, indici_bande(i,2)] = min(abs(f_9 - bande(i,2)));
end
indici_bande = min(max(indici_bande, 1), length(Pxx_9));
disp(indici_bande)
    Ptot_9=sum(Pxx_9);
    P_9=sum(Pxx_9(indici_bande(1,1):indici_bande(1,2))); 
    Pd_9=P_9/Ptot_9;
    P_9=sum(Pxx_9(indici_bande(2,1):indici_bande(2,2)));
    Pt_9=P_9/Ptot_9;
    P_9=sum(Pxx_9(indici_bande(3,1):indici_bande(3,2)));
    Pa_9=P_9/Ptot_9;
    P_9=sum(Pxx_9(indici_bande(4,1):indici_bande(4,2)));
    Pb1_9=P_9/Ptot_9;
    P_9=sum(Pxx_9(indici_bande(5,1):indici_bande(5,2)));
    Pb2_9=P_9/Ptot_9; 
    fprintf('\n Pd_9 = %.4f Pt_9 = %f Pa_9 = %f Pb1_9 = %f Pb2_9 = %f \n',Pd_9,Pt_9,Pa_9,Pb1_9,Pb2_9);


    %% correlogramma
    %partendo da 512 valori della sequenza di autocorrelazione finestrati mediante una finestra rettangolare, di 
    % Hamming, di Kaiser (512,5) e di Kaiser (512,12)
    N_corr=512;
    Q=1;

    %finestra rettangolare canale 7, Q=1
    w1=rectwin(N_corr);
    x=canale7-mean(canale7);
    [Pxxc,f] = correlogramma(N_corr,Q,w1,fc,x);
    figure(1)
    plot(f,abs(Pxxc ./ max(abs(Pxxc))), 'g');
    ris=input('inserire valore risoluzione: '); % se non data una ris_app da testo, è = a teorica cioè fc/N_corr
    title(['Correlogramma (Ris. ' num2str(ris) ' Hz)']);

    %finestra rettangolare canale 9, Q=1
    x=canale9-mean(canale9);
    [Pxxc,f] = correlogramma(N_corr,Q,w1,fc,x);
    figure(1)
    plot(f,abs(Pxxc ./ max(abs(Pxxc))), 'g');
    ris=input('inserire valore risoluzione: '); % se non data una ris_app da testo, è = a teorica cioè fc/N_corr
    title(['Correlogramma (Ris. ' num2str(ris) ' Hz)']);

    %finestra hamming canale 7, Q=1
    w1=hamming(N_corr/4);
    x=canale7-mean(canale7);
    [Pxxc,f] = correlogramma(N_corr,Q,w1,fc,x);
    figure(3)
    plot(f, abs(Pxxc ./ max(abs(Pxxc))), 'g');
    ris=input('inserire valore risoluzione: '); % se non data una ris_app da testo, è = a teorica cioè fc/N_corr
    title(['Correlogramma (Ris. ' num2str(ris) ' Hz)']);

    %finestra hamming Canale 9, Q=1
    x=canale9-mean(canale9);
    [Pxxc,f] = correlogramma(N_corr,Q,w1,fc,x);
    figure(4)
    plot(f, abs(Pxxc ./ max(abs(Pxxc))), 'g');
    ris=input('inserire valore risoluzione: '); % se non data una ris_app da testo, è = a teorica cioè fc/N_corr
    title(['Correlogramma (Ris. ' num2str(ris) ' Hz)']);

    %finestra Kaiser (512,5), canale 7, Q=1
    w1=kaiser(N_corr,5);
    x=seg7-mean(seg7);
    [Pxxc,f] = correlogramma(N_corr,Q,w1,fc,x);
    figure(5)
    plot(f, abs(Pxxc ./ max(abs(Pxxc))), 'g');
    ris=input('inserire valore risoluzione: '); % se non data una ris_app da testo, è = a teorica cioè fc/N_corr
    title(['Correlogramma (Ris. ' num2str(ris) ' Hz)']);

    %finestra Kaiser(512,5) Canale 9, Q=1
    x=seg9-mean(seg9);
    [Pxxc,f] = correlogramma(N_corr,Q,w1,fc,x);
    figure(6)
    plot(f, abs(Pxxc ./ max(abs(Pxxc))), 'g');
    ris=input('inserire valore risoluzione: '); % se non data una ris_app da testo, è = a teorica cioè fc/N_corr
    title(['Correlogramma (Ris. ' num2str(ris) ' Hz)']);

    %% Q=0.05, fai correlogramma per i due canali
    %non essendo specificato metodo utilizzo pwelch con hamming
    N_corr=512;
    Q=0.05;
    %finestra hamming canale 7, Q=0.05
    w1=hamming(N_corr/4);
    x=canale7-mean(canale7);
    [Pxxc,f] = correlogramma(N_corr,Q,w1,fc,x);
    figure(7)
    plot(f, abs(Pxxc ./ max(abs(Pxxc))), 'g');
    ris=input('inserire valore risoluzione: '); % se non data una ris_app da testo, è = a teorica cioè fc/N_corr
    title(['Correlogramma (Ris. ' num2str(ris) ' Hz)']);
    Pxxc=Pxxc_7;
    f=f_7;
    %finestra hamming Canale 9, Q=0.05
    x=canale9-mean(canale9);
    [Pxxc,f] = correlogramma(N_corr,Q,w1,fc,x);
    figure(8)
    plot(f, abs(Pxxc ./ max(abs(Pxxc))), 'g');
    ris=input('inserire valore risoluzione: '); % se non data una ris_app da testo, è = a teorica cioè fc/N_corr
    title(['Correlogramma (Ris. ' num2str(ris) ' Hz)']);
    Pxxc=Pxxc_9;
    f=f_9;
    %% suddivisione della potenza nelle bande con Q=0.05
    %canale 7
    bande = [0.5, 3.5; 3.5, 7; 7, 14; 14, 21; 21, 30];
indici_bande = zeros(size(bande));
for i = 1:size(bande,1)
    [~, indici_bande(i,1)] = min(abs(f_7 - bande(i,1)));
    [~, indici_bande(i,2)] = min(abs(f_7 - bande(i,2)));
end
indici_bande = min(max(indici_bande, 1), length(Pxxc_7));
disp(indici_bande)

 Ptot=sum(Pxxc_7);
    P=sum(Pxxc_7(indici_bande(1,1):indici_bande(1,2))); 
    Pd_7=P/Ptot;
    P=sum(Pxxc_7(indici_bande(2,1):indici_bande(2,2)));
    Pt_7=P/Ptot;
    P=sum(Pxxc_7(indici_bande(3,1):indici_bande(3,2)));
    Pa_7=P/Ptot;
    P=sum(Pxxc_7(indici_bande(4,1):indici_bande(4,2)));
    Pb1_7=P/Ptot;
    P=sum(Pxxc_7(indici_bande(5,1):indici_bande(5,2)));
    Pb2_7=P/Ptot; 
    fprintf('\n Pd_7 = %.4f Pt_7 = %f Pa_7 = %f Pb1_7 = %f Pb2_7 = %f \n',Pd_7,Pt_7,Pa_7,Pb1_7,Pb2_7);

    %canale 9
     bande = [0.5, 3.5; 3.5, 7; 7, 14; 14, 21; 21, 30];
indici_bande = zeros(size(bande));
for i = 1:size(bande,1)
    [~, indici_bande(i,1)] = min(abs(f_9 - bande(i,1)));
    [~, indici_bande(i,2)] = min(abs(f_9 - bande(i,2)));
end
indici_bande = min(max(indici_bande, 1), length(Pxxc_9));
disp(indici_bande)

 Ptot=sum(Pxxc_9);
    P=sum(Pxxc_9(indici_bande(1,1):indici_bande(1,2))); 
    Pd_9=P/Ptot;
    P=sum(Pxxc_9(indici_bande(2,1):indici_bande(2,2)));
    Pt_9=P/Ptot;
    P=sum(Pxxc_9(indici_bande(3,1):indici_bande(3,2)));
    Pa_9=P/Ptot;
    P=sum(Pxxc_9(indici_bande(4,1):indici_bande(4,2)));
    Pb1_9=P/Ptot;
    P=sum(Pxxc_9(indici_bande(5,1):indici_bande(5,2)));
    Pb2_9=P/Ptot; 
    fprintf('\n Pd_9 = %.4f Pt_9 = %f Pa_9 = %f Pb1_9 = %f Pb2_9 = %f \n',Pd_9,Pt_9,Pa_9,Pb1_9,Pb2_9);


    %% mappa codificata colore o mappa valore efficace
    % Costruzione della matrice
    Z=[std(eeg(1:5,:)); std(eeg(6:10,:)); std(eeg(11:15,:)); std(eeg(16:20,:)); std(eeg(21:25,:))]';

    % Interpolazione
    figure(9)
    mappa(Z,'Mappa distribuzione valore efficace segnale EEG');

    %% mappe della potenza associata ai 5 ritmi
    for i=1:5
    for j=1:5
        x2=eeg(:,j+(i-1)*5);
        [Zd(i,j),Zt(i,j),Za(i,j),Zb1(i,j),Zb2(i,j),Zu(i,j)]=rel_pot(x2);
    end
end
figure(10)
mappa(Zd,'Mappa ritmo delta');
figure(11)
mappa(Zt,'Mappa ritmo teta');
figure(12)
mappa(Za,'Mappa ritmo alfa');
figure(13)
mappa(Zb1,'Mappa ritmo beta 1');
figure(14)
mappa(Zb2,'Mappa ritmo beta 2');

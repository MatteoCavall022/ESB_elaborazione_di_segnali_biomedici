%xor deve essere lungo almeno quanto il numero di campioni per ogni epoca
%(N), se è più lungoo non importa abbiamo già previsto come gestirlo, se è
%più corto sapendo N da testo concateni a xor tanti zeri quanti
%N-length(xor)

S = load('signorm5.mat');
x = S.x;      % o il nome effettivo della variabile nel .mat
xor = S.xor;  % idem

 subplot(2,1,1)
 plot(xor(1:300),'r'); %Plottare tre epoche del segnale originale, dal 1� al 300� campione
 axis([0 300 -0.3 0.3]);
 title('Synthesized potential')

 %--------------------
 % Dal segnale di riferimento, determinare i range (in campioni)
 % entro i quali calcolare:
 % - il segnale picco-picco 
 % - il rumore
 %--------------------
 limit_signal=input('inserire valore fino a dove arriva il segnale %d');
 range_pp = 1:limit_signal; % a occhio dal grafico vedi che fino a dove hai segnale
 range_rumore = limit_signal:input('inserire valore fino a dove hai solo rumore'); % da 27 a circa 100 solo rumore
 
 signal_pp=max(xor(range_pp))-min(xor(range_pp));% xor perchè segnale pulito
 noise_stdev=std(x(range_rumore)); % x perchè segnale corrotto da rumore
 
 %--------------------
 % Calcolare il rapporto segnale-rumore
 %--------------------
 
 SNRinit =signal_pp/(4*noise_stdev);
 
 subplot(2,1,2) % numero di righe, colonne, posizione del grafico attivo
 plot(x(1:300),'g');
 axis([0 300 -0.3 0.3]);
 title(['Corrupted potential - SNR = ', num2str(SNRinit,'%4.1f')]) %num2str(SNRinit,'%4.1f' comando per fare in modo che a destra dell'=ti venga valore di SNRinit

 
 
%input per la funzione
N=100;
jitter=0;


 [SNR_teor, SNR_sper, xav] = averaging(x,xor, N, jitter, SNRinit, range_rumore, signal_pp);

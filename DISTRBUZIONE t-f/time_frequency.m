function [w,cw,aa,pxx] = time_frequency(x,ntlag,T,fc)
% x, segnale in input
% ntlag, numero di ritardi di solito = a N/2
% fc, frequenza campionamento
%T, supporto temporale o durata del segnale/brano

N=fc*T;
asset = (0:N-1)/fc;
xa=hilbert(x);
% Wigner-Ville
w = wig(xa, fc,ntlag); 
% Choi_williams
sigma=input('inserisci valore sigma: '); %puoi variare sigma per vedere un effetto filtrante differente
cw = dcw(xa, fc,ntlag, sigma,T); %calcola la act,funzione di ambiguità,la moltoiplica per kernel e ritrasforma sul piano t-f
% commento sui valori di sigma:  Sigma piccolo: più filtraggio delle interferenze, ma risoluzione tempo-frequenza minore.
% Sigma grande: meno filtraggio, risoluzione maggiore, ma più interferenze
%% plot del segnale 
% Wigner-Ville
upper = fc/4;
[roww, ~] = size(w);
assef = (0:roww-1)/roww*upper;

figure (1);
mesh(asset, assef, real(w));
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Wigner-Ville: Segnale ');
axis([min(asset) max(asset) min(assef) max(assef) 2*min(min(real(w))) 2*max(max(real(w)))]);
figure(2);
contour(asset, assef, real(w),20);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Wigner-Ville: Segnale ');
axis([min(asset) max(asset) min(assef) max(assef) 2*min(min(real(w))) 2*max(max(real(w)))]);
% Choi-williams
%sigma 1
figure (3);
mesh(asset, assef, real(cw));
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title(['Choi-Williams: Segnale  (\sigma= ' num2str(sigma) ')']);
axis([min(asset) max(asset) min(assef) max(assef) 2*min(min(real(cw))) 2*max(max(real(cw)))]);
figure(4);
contour(asset, assef, real(cw),20);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Choi-Williams: Segnale ');
axis([min(asset) max(asset) min(assef) max(assef) 2*min(min(real(cw))) 2*max(max(real(cw)))]);
%% calcolo funzione ambiguità
aa = daf(xa, ntlag, fc, T); % xa1 = segnale analitico, ntlag, fc, T

% aa1 è una matrice complessa: modulus = ampiezza della funzione di ambiguità
figure;
imagesc(abs(aa));
xlabel('Delay (samples)');
ylabel('Doppler (samples)');
title('Discrete Ambiguity Function (DAF) - Segnale 1');
colorbar;
%% Marginale in frequenza: somma sulle colonne (asse tempo) e periodogramma
marginale_f = sum(abs(w), 2);
upper=fc/2;
assef = (0:roww-1)/roww*upper;
% Periodogramma semplice
NFFT=N; % altrimenti NFFT=fc/ris_app
[pxx, f] = pwelch(x-mean(x), rectwin(N), 0, NFFT,fc);

figure;
plot(assef, marginale_f / max(marginale_f), 'b', 'LineWidth', 1.5); hold on;
plot(f, pxx / max(pxx), 'r--', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Normalized amplitude');
legend('Marginale in frequenza', 'Periodogramma');
title('Marginale vs Periodogramma - Segnale');
end
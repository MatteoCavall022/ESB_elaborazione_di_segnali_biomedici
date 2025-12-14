function [Pxxc,f] = correlogramma(N_corr,Q,w1,fc,x)
%N_corr è dato da testo, valori della sequenza di autocorrelazione, se non
%dato da testo ris=1/((N_corr/fc)*Q);%N_corr/fc=T
%Q, se non dice niente metti =1
%w1 è il tipo di finestra (rectwin, hamming, keiser) e la scrivi 
% w1=rectwin[N_corr], w1=kaiser(N_corr,5)
% x è pezzo di segnale senza valor medio, se non fatto prima FAI SUBITO
% X-MEAN(X) prima di richiamare la funzione
%calcolo correlogramma
disp('Hai tolto il valore medio ad X?')
if Q==1
    risapp=input('inserire valore risapp da testo');
    NFFT=floor(fc/risapp);
else 
    NFFT=floor(N_corr*Q);
end
% se Q=1 e ti dice che vuole un particolare valore di risapp devi usare quella per calcolare NFFT, altrimenti, e se varia Q (in tal caso sempre!), 
% risapp=teorica e quindi NFFT=N_corr*Q cioè quello che trovi per la ris teorica*Q oppure NFFT=fc/risteorica  
ntlag = round((N_corr-1)/2); 
acs = xcorr(x, ntlag, 'biased');
acs1 = ntlag + 1 - (length(w1)/2) : ntlag + (length(w1)/2);
acsw = acs(acs1) .* w1;
Pxxc = fft(acsw, NFFT); %gira acsw
%elimina la replica spettrale, cioè prendi la primma metà
Pxxc = Pxxc(1:(length(Pxxc)/2)+1);
f = fc * (0:(length(Pxxc)-1)) / NFFT;
end
%fuori dalla funzione fai plot(f, abs(Pxxc ./ max(abs(Pxxc))), 'g');
% title(['Correlogramma (Ris. ' num2str(ris) ' Hz)']);

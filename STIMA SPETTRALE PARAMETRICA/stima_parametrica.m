function [LF_HF_ratio,a] = stima_parametrica(xm,NFFT,ord_iniziale, ord_finale,metodo)
% fc, frequenza campionamento
%xm, segnale RR 
%NFFT, da dare in ingresso a pburg, se ti da ris_teorica trovi
%NFFT=fc/ris_teorica
%ord_iniziale e ord_finale stabiliscono intervallo di ordini per calcolare
%la varianza dell'errore
%metodo, puoi scegliere burg,covarianza_mod,yule_walker

aumento=0.05; %valore in % di cui devo aumentare la varianza dell'errore per trovare il minimo rispetto 
% alla var asintotica. nel pdf era Scegliere come ordine del modello il minimo valore per cui
%la varianza dell'errore risulta non superiore alla varianza asintotica  aumentata del 5%.
xm2=xm;
xm=xm-mean(xm);

if strcmp(metodo, 'burg')
%calcolo varianza errore
for NN = ord_iniziale:ord_finale
    [a,e(NN)] = arburg(xm,NN); %devo dare in ingresso a aburg segnale e ordine e in uscita da vettore dei coefficienti a stimati e la varianza dell'errore di predizione
end
% NB: la varianza asintotica sarà l'ultimo valore del vettore contenente la
% varianza + aumento
asint =(e(ord_finale).*aumento)+ e(ord_finale);

% Trovare gli indici dove ci sono valori di varianza che sono minori della varianza asintotica

ind= find(e < asint);    


% Scegliere l'ordine; 
%Scegliere come ordine del modello il minimo valore per cui
%la varianza dell'errore risulta non superiore alla varianza asintotica 
% Si consiglia di considerare il secondo elemento perché il primo contiene un valore inesatto
order =ind(2);

% Stima spettrale parametrica con metodo di Burg.
%[pxx, f] = pburg(x, p, nfft); con p ordine del modello AR
[Pb, w] = pburg(xm,order,NFFT);

disp('il segnale RR era in ms e volevo output in s per cui ho diviso per 1000 asse f se hai già in s togli 1000')
f = w/(2*pi)/(mean(xm2)/1000); %normalizzazione delle frequenze
figure, plot(f, Pb/max(Pb)), title('PSD di Burg'), xlabel('frequenza (cicli/s)')
axis tight


%calcolare il rapporto LF/HF dividendo la potenza dello spettro AR nella
% banda LF (60 - 140 mHz) per quella nella banda HF (140 - 300 mHz). Tenere conto del
% fatto che gli estremi della banda HF potrebbero variare a seconda dello stato del soggetto
% e, quindi, del ritmo respiratorio.
LF_min = 0.060;   % 60 mHz
LF_max = 0.140;   % 140 mHz
HF_min = 0.140;   % 140 mHz
HF_max = 0.300;   % 300 mHz

% Trova gli indici delle bande
indLF = find(f >= LF_min & f < LF_max);
indHF = find(f >= HF_min & f <= HF_max);

% Calcola la potenza nelle due bande come somma della PSD (area)
P_LF = sum(Pb(indLF));
P_HF = sum(Pb(indHF));

% Calcola il rapporto LF/HF
LF_HF_ratio = P_LF / P_HF;

elseif strcmp(metodo, 'covarianza_mod')
    for NN = ord_iniziale:ord_finale
    [a,e(NN)] = armcov(xm,NN); %devo dare in ingresso a armcov segnale e ordine e in uscita da vettore dei coefficienti a stimati e la varianza dell'errore di predizione
end

% Calcolo della varianza asintotica per ordini del modello tra 2 e 50
% Calcolare il 5% della varianza asintotica e sommare alla varianza
% asintotica
asint =(e(ord_finale).*aumento)+ e(ord_finale);
    
% Trovare gli indici in cui i valori di varianza che sono minori della varianza asintotica
% aumentata del 5%
ind=find(e<asint);
    
% Scegliere l'ordine; 
% Si consiglia di considerare il secondo elemento perché il primo contiene un valore inesatto
order=ind(2);
    
% Stima spettrale parametrica con metodo della Covarianza Modificata.
%[pxx, f] = pmcov(x, p, nfft); con p ordine del modello AR
[Pb, w] = pmcov(xm,order,NFFT);
% Normalizzare le frequenze
% Rappresentare graficamente la densità spettrale di potenza della serie RR
% per l'ordine ottenuto con il metodo della covarianza modificata
disp('il segnale RR era in ms e volevo output in s per cui ho diviso per 1000 asse f se hai già in s togli 1000')
f = w/(2*pi)/(mean(xm2)/1000); %normalizzazione delle frequenze
figure, plot(f, Pb/max(Pb)), title('PSD di Covarianza modificata'), xlabel('frequenza (cicli/s)')
axis tight
  
LF_min = 0.060;   % 60 mHz
LF_max = 0.140;   % 140 mHz
HF_min = 0.140;   % 140 mHz
HF_max = 0.300;   % 300 mHz

% Trova gli indici delle bande
indLF = find(f >= LF_min & f < LF_max);
indHF = find(f >= HF_min & f <= HF_max);

% Calcola la potenza nelle due bande come somma della PSD (area)
P_LF = sum(Pb(indLF));
P_HF = sum(Pb(indHF));

% Calcola il rapporto LF/HF
LF_HF_ratio = P_LF / P_HF;

elseif strcmp(metodo, 'yule_walker')

    for NN = ord_iniziale:ord_finale
    [a,e(NN)] = aryule(xm,NN); %devo dare in ingresso a armcov segnale e ordine e in uscita da vettore dei coefficienti a stimati e la varianza dell'errore di predizione
end

% Calcolo della varianza asintotica per ordini del modello tra 2 e 50


% Calcolare l'aumento della varianza asintotica e sommare alla varianza
% asintotica
asint =(e(ord_finale).*aumento)+ e(ord_finale);
    
% Trovare gli indici dei valori di varianza che sono minori della varianza asintotica
% aumentata di aumento
ind=find(e<asint);
    
% Scegliere l'ordine; 
% Si consiglia di considerare il secondo elemento perché il primo contiene 
% un valore inesatto
order=ind(2);
    
% Stima spettrale parametrica con metodo yulewalker
%[pxx, f] = pmcov(x, p, nfft); con p ordine del modello AR
[Pb, w] = pyulear(xm,order,NFFT);
% Normalizzare le frequenze
% Rappresentare graficamente la densità spettrale di potenza della serie RR
% per l'ordine ottenuto con il metodo yulewalker
disp('il segnale RR era in ms e volevo output in s per cui ho diviso per 1000 asse f se hai già in s togli 1000')
f = w/(2*pi)/(mean(xm2)/1000); %normalizzazione delle frequenze, per quwllo sopra dopo NFFT non metti fc perchè qui hai i battiti non una fcapionamento 
figure, plot(f, Pb/max(Pb)), title('PSD di Yulewalker'), xlabel('frequenza (cicli/s)')
axis tight

LF_min = 0.060;   % 60 mHz
LF_max = 0.140;   % 140 mHz
HF_min = 0.140;   % 140 mHz
HF_max = 0.300;   % 300 mHz

% Trova gli indici delle bande
indLF = find(f >= LF_min & f < LF_max);
indHF = find(f >= HF_min & f <= HF_max);

% Calcola la potenza nelle due bande come somma della PSD (area)
P_LF = sum(Pb(indLF));
P_HF = sum(Pb(indHF));

% Calcola il rapporto LF/HF
LF_HF_ratio = P_LF / P_HF;

end
end
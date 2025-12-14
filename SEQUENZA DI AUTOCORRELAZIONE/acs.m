function xacs = acs(x,ntlag,bias)
%Function xacs = acs(x,ntlag,bias);
%
% Funzione per la stima della sequenza di autocorrelazione (ACS)
% di una sequenza numerica di ingresso 'x'. La variabile 'ntlag'
% determina i ritardi temporali su cui stimare la ACS. La lunghezza
% della sequenza di uscita ? pari a 2*ntlag+1.
% Il numero massimo di ritardi temporali possibili per la stima
% ? pari alla lunghezza di 'x'.
%
% La funzione ? in grado di effettuare la stima della ACS sia mediante
% lo stimatore polarizzato che mediante quello non polarizzato. La
% variabile logica 'bias' determina quale stimatore utilizzare.
%
% Ingressi:     x       vettore contenente la sequenza di cui stimare l'ACS
%               ntlag   numero di ritardi temporali per la stima
%               bias    variabile logica per la scelta dello stimatore
%                       se bias = 0 -> stimatore non polarizzato
%                       se bias = 1 -> stimatore polarizzato

disp('in input alla funzione devi dare x togliendo prima il valor medio')

%ciclo principale che conta i ritardi
for i = 0:ntlag
    if i == 0
        %stima ACS per ritardo nullo
        xacs(i+1)= sum(x.* conj(x)); %i+1 è m, se i=0 (ritardo nullo) entra in questa riga e con .* moltiplica el per el tutti i valori del vettore, ottengo un vettore e faccio sum perchè nella formula voglio sommatoria
    % quando ritardo è nullo else sotto non lo vede quindi devo
    % normalizzare subito 
    xacs=xacs/length(x);
    else
        %stima ACS per ritardi positivi
        xacs(i+1)=sum([zeros(1,i), x(1:end-i)].*conj(x));

        %stima ACS per ritardi negativi, qui mette (se ntlag=5 per es) il
        %valore sarà in posizione 10 se invece vale 1 in posizione 3 quindi
        %da dx verso sx poi con fftshift lo gira da sx a dx
        xacs(2*ntlag-i+2)= conj(sum([zeros(1,i), x(1:end-i)].*conj(x))); %sfrutto proprietà che rxx[-m]=compl con rxx[m]
        
        %normalizzazione della ACS stimata
        if bias == 0
            xacs(i+1) = xacs(i+1)/(length(x)-i);
            xacs(2*ntlag-i+2) = xacs(2*ntlag-i+2)/(length(x)-i); %avevamo messo abs i come da formula però in realtà qui divide già in base a ritardo positivo e negativo quindi non serve
        else
            xacs(i+1) = xacs(i+1)/length(x);
            xacs(2*ntlag-i+2) =xacs(2*ntlag-i+2)/length(x);
        end
    end
end

xacs = fftshift(xacs);    
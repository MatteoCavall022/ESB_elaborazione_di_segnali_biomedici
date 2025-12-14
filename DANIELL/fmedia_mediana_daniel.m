function [fmean, fmed] = fmedia_mediana_daniel(Pd, f)
% Calcola la frequenza media e mediana da uno spettro levigato (Daniell)
% Pd = spettro di potenza
% f = vettore frequenze (stesso ordine di Pd)

% Rimozione eventuali NaN (per sicurezza)
valid = ~isnan(Pd) & ~isnan(f);
Pd = Pd(valid);
f = f(valid);

% Frequenza media
fmean = sum(Pd .* f) / sum(Pd);

% Frequenza mediana
df = f(2) - f(1);  % passo frequenziale costante
area_tot = sum(Pd) * df;
area_cum = cumsum(Pd) * df;
idx = find(area_cum >= area_tot / 2, 1, 'first');
fmed = f(idx);
end

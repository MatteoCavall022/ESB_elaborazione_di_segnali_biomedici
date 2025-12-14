function [rms_vals, fmean_vals, fmedian_vals, cv_vals] = fatigue_plot(sig, fs, epoch_len_sec, d)
% FATIGUE_PLOT - Calcola e grafica i fatigue plot di un segnale EMG multicanale
% in funzione di RMS, MNF, MDF e velocità di conduzione (VC)
%
% INPUT:
%   sig            = matrice EMG (canali x campioni), es. [3 x N]
%   fs             = frequenza di campionamento [Hz]
%   epoch_len_sec  = durata epoche in secondi (es. 0.5)
%   d              = distanza inter-elettrodica (m), es. 0.005 m (5 mm)

[num_channels, N] = size(sig);
if num_channels < 3
    error('Sono necessari almeno 3 canali EMG (righe)');
end

hop = round((epoch_len_sec * fs) / 2);   % 50% overlap
epoch_len = round(epoch_len_sec * fs);
starts = 1:hop:(N - epoch_len + 1);

% Output
rms_vals = zeros(1, length(starts));
fmean_vals = zeros(1, length(starts));
fmedian_vals = zeros(1, length(starts));
cv_vals = zeros(1, length(starts));
time_vals = (starts + epoch_len/2) / fs;

for i = 1:length(starts)
    idx = starts(i):(starts(i) + epoch_len - 1);
    epoch_dd1 = sig(1, idx) - sig(2, idx);  % Canale 1 - Canale 2
epoch_dd2 = sig(2, idx) - sig(3, idx);  % Canale 2 - Canale 3

    epoch_central = sig(2, idx);

    % RMS
    rms_vals(i) = sqrt(mean(epoch_central.^2));

    % Fmean e Fmedian con funzione esterna
    [fmean_vals(i), fmedian_vals(i)] = fmean(epoch_central, fs, epoch_len_sec);

    % CV con metodo spettrale (via delay)
    fft_dd1 = fft(epoch_dd1 - mean(epoch_dd1), epoch_len);
    fft_dd2 = fft(epoch_dd2 - mean(epoch_dd2), epoch_len);

    fft1r = real(fft_dd1);
    fft1i = imag(fft_dd1);
    fft2r = real(fft_dd2);
    fft2i = imag(fft_dd2);

    try
        delay_samples = delay(fft1r, fft1i, fft2r, fft2i, 0);
        delay_time = delay_samples / fs;
        if abs(delay_time) > 0
            cv_vals(i) = d / abs(delay_time);
        else
            cv_vals(i) = NaN;
        end
    catch
        cv_vals(i) = NaN;
    end
end

% Plot
figure;
subplot(4,1,1);
plot(time_vals, rms_vals, 'b-o');
ylabel('RMS [mV]'); title('RMS nel tempo'); grid on;

subplot(4,1,2);
plot(time_vals, fmean_vals, 'r-o');
ylabel('MNF [Hz]'); title('Frequenza media nel tempo'); grid on;

subplot(4,1,3);
plot(time_vals, fmedian_vals, 'm-o');
ylabel('MDF [Hz]'); title('Frequenza mediana nel tempo'); grid on;

subplot(4,1,4);
plot(time_vals, cv_vals, 'g-o');
ylabel('VC [m/s]'); xlabel('Tempo [s]'); title('Velocità di conduzione'); grid on;

sgtitle('Fatigue Plot EMG - RMS, MNF, MDF, CV');
end

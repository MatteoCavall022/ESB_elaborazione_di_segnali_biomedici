function [Pd, f] = daniel(x, P, fs)
    NFFT = 5 * fs;  % segnale lungo 5s → ris = 0.2 Hz
    [Pxx, f] = pwelch(x - mean(x), boxcar(NFFT), 0, NFFT, fs);
    Pd = zeros(size(Pxx));
    for i = 1:length(Pxx)
        if i <= P
            Pd(i) = mean(Pxx(1:P+i));
        elseif i > length(Pxx)-P
            Pd(i) = mean(Pxx(i-P:end));
        else
            Pd(i) = mean(Pxx(i-P:i+P));
        end
    end
end



function [SNR_teor, SNR_sper, xav] = averaging_corrected(x,xor, N, jitter, SNRinit, range_rumore, signal_pp)
%xor: segnale pulito lungo 1 epoca, vettore colonna
%x: segnale sporco completo, vettore colonna
%N: numero campioni per epoca
%jitter: se presente =1, altrimenti 0

n_epoch = floor(length(x)/N);
disp ('se hai jitter hai messo in cartella la funzione delay?')

if jitter == 1
    xor = xor(1:N);
    for i = 0:(n_epoch - 2)
        xact = x((N+1)+N*i:(2*N)+N*i);
        xc = xcorr(xact, xor);  % cross-correlazione
        M = max(xc);
        indM = find(xc == M) - N;
        d(i+2) = delay(real(fft(xact)), imag(fft(xact)), real(fft(xor)), imag(fft(xor)), indM);
    end
    d = round(d);

    xrial(1,:) = x(1:N);
    for i = 0:(n_epoch - 2)
        segment = x(i*N + (N+1):i*N + (2*N));
        if d(i+2) == 0
            xrial(i+2,:) = segment;
        else
            xrial(i+2,:) = circshift(segment, -d(i+2));
        end
    end
end

n_epoch = fix(length(x)/N);
FR = n_epoch;

if jitter == 1
    x_matrix = xrial';
    FR = size(x_matrix, 2);
else
    x_matrix = reshape(x, N, n_epoch);
end

sum_epoch = zeros(N,1);
SNR_teor = zeros(FR,1);
SNR_sper = zeros(FR,1);

for jj = 1:FR
    sum_epoch = sum_epoch + x_matrix(:,jj);
    xav_jj = sum_epoch / jj;

    noise_stdev = std(xav_jj(range_rumore));
    SNR_teor(jj) = 20 * log10(SNRinit * sqrt(jj));
    SNR_sper(jj) = 20 * log10(signal_pp / (4 * noise_stdev));

    subplot(2,1,1)
    plot(xor(1:N), 'r')
    axis([0 N -0.3 0.3])
    title('Synthesized potential')

    subplot(2,1,2)
    plot(xav_jj(1:100), 'g')
    axis([0 N -0.3 0.3])
    title(['Averaged potential - Epoch ', int2str(jj), '    SNR = ', num2str(SNR_sper(jj), '%4.1f')])
    drawnow

    t1 = clock;
    while etime(clock, t1) <= 0.5
    end
end

xav = sum_epoch / FR;

figure(3)
plot(SNR_teor, 'r')
hold on
plot(SNR_sper, 'b')
xlabel('n_epoch')
ylabel('SNR (in dB)')
end

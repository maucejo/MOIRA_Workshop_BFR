function y = agwn(x, SNR_dB)
    %y = awgn_noise(x,SNR) adds AWGN noise vector to signal 'x' to generate a
    %resulting signal vector y of specified SNR in dB

    % rng('default');               %set the random generator seed to default (for comparison only)
    [N, L] = size(x);
    SNR = 10^(SNR_dB/10);         % SNR to linear scale
    Esym = sum(abs(x).^2, 2)/(L); % Calculate actual symbol energy
    N0 = Esym/SNR;                % Find the noise spectral density
    if(isreal(x))
        noiseSigma = repmat(sqrt(N0), 1 ,L);      %Standard deviation for AWGN Noise when x is real
        n = noiseSigma.*randn(N, L); %computed noise
    else
        noiseSigma = repmat(sqrt(N0/2), 1, L);                       %Standard deviation for AWGN Noise when x is complex
        n = noiseSigma.*(randn(N, L) + 1i*randn(N, L)); %computed noise
    end
    y = x + n; %received signal
end

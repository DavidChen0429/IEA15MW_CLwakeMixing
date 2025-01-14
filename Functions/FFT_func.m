function [f, P1] = FFT_func(signal, fil_index, Fs)
    filtered_signal = signal(fil_index:end);   % filter out nonexcited period
    filtered_signal = filtered_signal - mean(filtered_signal);
    T = 1/Fs;
    L = length(filtered_signal);
    n = 2^nextpow2(L);         % Next prime2 number larger than L for faster
    f = Fs/n * (0:(n/2));
    signal_FFT = fft(filtered_signal, n);
    P2 = abs(signal_FFT/n);    % normalization
    P1 = P2(1:n/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

%     figure()
%     plot(f,P1,"LineWidth",1);
%     title("Filtered signal in Frequency Domain")
%     xlabel("f (Hz)")
%     ylabel("Magnitude")
%     text_position = [0.95, 0.05]; % Normalized coordinates (bottom right)
%     annotation_text = sprintf('FFT information');
%     text('Units', 'normalized', 'Position', text_position, ...
%         'String', annotation_text, 'HorizontalAlignment', ...
%         'right', 'VerticalAlignment', 'bottom');
end
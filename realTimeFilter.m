% Example usage
Fs = 1000; % Sampling frequency
Fc = 100;  % Cutoff frequency

% Define the filter parameters
Wn = Fc / (Fs / 2);
[b, a] = butter(4, Wn, 'low');

% Compute the frequency response
[h, w] = freqz(b, a, 1024, Fs);

% Plot the magnitude response
figure;
subplot(2, 1, 1);
plot(w, 20*log10(abs(h)));
title('Butterworth Filter Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

% Plot the phase response
subplot(2, 1, 2);
plot(w, angle(h) * (180/pi));
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
grid on;
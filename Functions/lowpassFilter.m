function [Signal_filtered] = lowpassFilter(Signal, Fs, Fc)
Wn = Fc / (Fs / 2);
[b, a] = butter(4, Wn, 'low');
Signal_filtered = filtfilt(b, a, Signal);
end


function [b_fir, n] = FIR_LPF(Fs,Fc)
Wn = Fc / (Fs / 2);
n = 50; % Filter order
b_fir = fir1(n, Wn, 'low');
end
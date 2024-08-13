N = 300;
AMPL = 1;
Ts = 1;
F = 0.5;
Fstop = 10;
a = idprbs(N, AMPL, Ts, F, Fstop);
plot(a)
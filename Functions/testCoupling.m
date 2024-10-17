function [RGA_ss, RGA_bw] = testCoupling(sys, bw, timeStep)
% bw is given in [Hz]
G_ss= dcgain(sys);
RGA_ss = G_ss .* (inv(G_ss))';

[G_mag, G_phase] = bode(sys, bw*2*pi);  
G_bw = evalfr(sys, exp(j*bw*2*pi*timeStep));
RGA_bw = abs(G_bw) .* (inv(abs(G_bw)))';
disp('==================================')
disp('RGA of steady-state:')
disp(RGA_ss)    % steady-state

disp('RGA of bandwidth frequency:')
disp(RGA_bw)    % bandwidth frequency
end
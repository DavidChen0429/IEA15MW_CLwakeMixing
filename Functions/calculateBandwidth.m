function [bandwidth_frequency] = calculateBandwidth(sys)
% Return the bandwidth in rad/s, if converting to Hz, don't forget to
% divide the value by 2*pi !!!!!!!!!!!!!!
    [mag, phase, w] = bode(sys);
    mag_dB = 20*log10(squeeze(mag));
    max_mag_dB = max(mag_dB);
    idx_bandwidth = find(mag_dB <= max_mag_dB-3, 1, 'first');
    bandwidth_frequency = w(idx_bandwidth);
    bandwidth_frequency = bandwidth_frequency / (2*pi); % Conver to Hz
    % fprintf('System bandwidth: %.5f Hz\n', bandwidth_frequency);
end
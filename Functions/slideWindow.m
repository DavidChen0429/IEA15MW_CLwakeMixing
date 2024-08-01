function [meanArray] = slideWindow(data,windowLength,delay)
Fs = 10;
Fc = 0.02;
buffer = [];
meanArray = [];
for i = 1:1:length(data)
    if i > delay
        buffer = [buffer data(i)];
        if length(buffer) == windowLength
            buffer = lowpassFilter(buffer, Fs, Fc);
            meanArray = [meanArray mean(buffer)];
            buffer = [];   % reset the sliding window
        end
    end
end
end
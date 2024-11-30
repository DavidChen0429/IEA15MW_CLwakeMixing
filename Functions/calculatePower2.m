function [PowerAvg] = calculatePower2(filter,Power_store)
    PowerAvg1 = mean(Power_store(filter:filter+3000, :))/1e3;
    PowerAvg2 = mean(Power_store(filter+1000:filter+4000, :))/1e3;
    PowerAvg3 = mean(Power_store(filter+2000:filter+5000, :))/1e3;
    PowerAvg4 = mean(Power_store(filter+3000:filter+6000, :))/1e3;
    PowerAvg = [PowerAvg1 PowerAvg2 PowerAvg3 PowerAvg4];
end
function [meanU, TI] = calculateTI(Datablock)
    windspeedData = Datablock.LiDAR_data;
    meanU = arrayfun(@(x) mean(x.u_los), windspeedData);

    TI = zeros(length(windspeedData), 1);
    for i = 1:length(windspeedData)
        TI(i) = std(windspeedData(i).u_los)/meanU(i);
    end
end
function [DEL] = calculateDEL(filter,Moop1_store, timeStep)
    Moop1_store_buf = Moop1_store(filter:end, :);
    cycles = rainflow(Moop1_store_buf, 1/timeStep);
    neq = 1;
    inverseWohlerSlope = 10;
    DEL = (sum(cycles(:, 1) .* (cycles(:, 2).^inverseWohlerSlope))/neq)^(1/inverseWohlerSlope);
end
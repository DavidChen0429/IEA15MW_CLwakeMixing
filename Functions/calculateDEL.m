function [DEL] = calculateDEL(Moop1_store, timeStep)
    cycles = rainflow(Moop1_store, timeStep);
    neq = 1;
    inverseWohlerSlope = 10;
    DEL = (sum(cycles(:, 1) .* (cycles(:, 2).^inverseWohlerSlope))/neq)^(1/inverseWohlerSlope);
end
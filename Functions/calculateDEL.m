function [DEL] = calculateDEL(filter,Mflap1_store,Medge1_store,Mflap2_store,Medge2_store,Mflap3_store,Medge3_store, timeStep)
    % This function computes the DEL of blades (edgewise and flapwise)
    neq = 1;
    inverseWohlerSlope = 10;
    fs = 1/timeStep;

    % Blade1
        % 1. flapwise 
    Mflap1_store_buf = Mflap1_store(filter:end, :);
    cycles_f1 = rainflow(Mflap1_store_buf, fs);
    DEL_f1 = (sum(cycles_f1(:, 1) .* (cycles_f1(:, 2).^inverseWohlerSlope))/neq)^(1/inverseWohlerSlope);
    DEL_f1 = (cycles_f1(:,2)'.^inverseWohlerSlope*cycles_f1(:,1)/neq)^(1/inverseWohlerSlope);
        % 2. edgewise
    Medge1_store_buf = Medge1_store(filter:end, :);
    cycles_e1 = rainflow(Medge1_store_buf, fs);
    DEL_e1 = (sum(cycles_e1(:, 1) .* (cycles_e1(:, 2).^inverseWohlerSlope))/neq)^(1/inverseWohlerSlope);
    DEL_e1 = (cycles_e1(:,2)'.^inverseWohlerSlope*cycles_e1(:,1)/neq)^(1/inverseWohlerSlope);
    
    % Blade2
        % 1. flapwise 
    Mflap2_store_buf = Mflap2_store(filter:end, :);
    cycles_f2 = rainflow(Mflap2_store_buf, fs);
    DEL_f2 = (sum(cycles_f2(:, 1) .* (cycles_f2(:, 2).^inverseWohlerSlope))/neq)^(1/inverseWohlerSlope);
    DEL_f1 = (cycles_f2(:,2)'.^inverseWohlerSlope*cycles_f2(:,1)/neq)^(1/inverseWohlerSlope);
        % 2. edgewise
    Medge2_store_buf = Medge2_store(filter:end, :);
    cycles_e2 = rainflow(Medge2_store_buf, fs);
    DEL_e2 = (sum(cycles_e2(:, 1) .* (cycles_e2(:, 2).^inverseWohlerSlope))/neq)^(1/inverseWohlerSlope);
    DEL_e2 = (cycles_e2(:,2)'.^inverseWohlerSlope*cycles_e2(:,1)/neq)^(1/inverseWohlerSlope);

    % Blade3
        % 1. flapwise 
    Mflap3_store_buf = Mflap3_store(filter:end, :);
    cycles_f3 = rainflow(Mflap3_store_buf, fs);
    DEL_f3 = (sum(cycles_f3(:, 1) .* (cycles_f3(:, 2).^inverseWohlerSlope))/neq)^(1/inverseWohlerSlope);
    DEL_f3 = (cycles_f3(:,2)'.^inverseWohlerSlope*cycles_f3(:,1)/neq)^(1/inverseWohlerSlope);
        % 2. edgewise
    Medge3_store_buf = Medge3_store(filter:end, :);
    cycles_e3 = rainflow(Medge3_store_buf, fs);
    DEL_e3 = (sum(cycles_e3(:, 1) .* (cycles_e3(:, 2).^inverseWohlerSlope))/neq)^(1/inverseWohlerSlope);
    DEL_e3 = (cycles_e3(:,2)'.^inverseWohlerSlope*cycles_e3(:,1)/neq)^(1/inverseWohlerSlope);

    % Array storage 
    DEL = [DEL_f1, DEL_e1, DEL_f2, DEL_e2, DEL_f3, DEL_e3];
    
end
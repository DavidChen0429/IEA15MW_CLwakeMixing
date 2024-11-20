function [PowerAvg] = calculatePower(filter,Power_store,D_NREL5MW,U_inflow)
%     Cp_store_buf = Cp_store(filter:end, :);
%     Power2 = (1.225*pi*((D_NREL5MW)/2)^2*U_inflow^3*Cp_store_buf)/2; % same Power_store
%     PowerAvg = mean(Power2)/1e6; % MW

    PowerAvg = mean(Power_store(filter:end, :))/1e3;
end
function [PowerAvg] = calculatePower(Cp_store,D_NREL5MW,U_inflow)
    Power2 = (1.225*pi*((D_NREL5MW)/2)^2*U_inflow^3*Cp_store)/2; % same Power_store
    PowerAvg = mean(Power2)/1e6; % MW
end
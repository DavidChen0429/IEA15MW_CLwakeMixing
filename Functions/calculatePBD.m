function [PBD] = calculatePBD(filter,PitchAngles,Mflap1_store,Medge1_store,Mflap2_store,Medge2_store,Mflap3_store,Medge3_store)
    % This function computes the PBD
    PitchAngles_buf = PitchAngles(filter:end, :);
    N = length(PitchAngles_buf);
    WohlerSlope = 3;
    gamma = 0:10:350;
    ThetaDiff = abs([PitchAngles_buf(1, :); diff(PitchAngles_buf)]);
    
    % Blade1
    PBD1 = 0;
    Mflap1_store_buf = Mflap1_store(filter:end, :);
    Medge1_store_buf = Medge1_store(filter:end, :);
    for k = 1:N
        buff1 = cosd(gamma)*Mflap1_store_buf(k)+sind(gamma)*Medge1_store_buf(k);
        max_buff1 = max(buff1);
        max_bufff1 = max(max_buff1, 0);
        PBD1 = PBD1 + ThetaDiff(k, 1)*max_bufff1^WohlerSlope;
    end
    PitchBearingDamage1 = PBD1 / 1e3; % [kNm deg]
    PitchBearningDEL1 = (PitchBearingDamage1 ./1e8).^(1/WohlerSlope);

    % Blade2
    PBD2 = 0;
    Mflap2_store_buf = Mflap2_store(filter:end, :);
    Medge2_store_buf = Medge2_store(filter:end, :);
    for k = 1:N
        buff2 = cosd(gamma)*Mflap2_store_buf(k)+sind(gamma)*Medge2_store_buf(k);
        max_buff2 = max(buff2);
        max_bufff2 = max(max_buff2, 0);
        PBD2 = PBD2 + ThetaDiff(k, 2)*max_bufff2^WohlerSlope;
    end
    PitchBearingDamage2 = PBD2 / 1e3; % [kNm deg]
    PitchBearningDEL2 = (PitchBearingDamage2 ./1e8).^(1/WohlerSlope);

    % Blade3
    PBD3 = 0;
    Mflap3_store_buf = Mflap3_store(filter:end, :);
    Medge3_store_buf = Medge3_store(filter:end, :);
    for k = 1:N
        buff3 = cosd(gamma)*Mflap3_store_buf(k)+sind(gamma)*Medge3_store_buf(k);
        max_buff3 = max(buff3);
        max_bufff3 = max(max_buff3, 0);
        PBD3 = PBD3 + ThetaDiff(k, 3)*max_bufff3^WohlerSlope;
    end
    PitchBearingDamage3 = PBD3 / 1e3; % [kNm deg]
    PitchBearningDEL3 = (PitchBearingDamage3 ./1e8).^(1/WohlerSlope);
    
    % Array storage
    PBD = [PitchBearningDEL1, PitchBearningDEL2, PitchBearningDEL3];
    
end
function [PitchBearingDamage] = calculatePBD(filter,PitchAngles,Mflap1_store,Medge1_store)
    PitchAngles_buf = PitchAngles(filter:end, :);
    Mflap1_store_buf = Mflap1_store(filter:end, :);
    Medge1_store_buf = Medge1_store(filter:end, :);
    N = length(PitchAngles_buf);
    WohlerSlope = 3;
    PBD = 0;
    gamma = 0:10:350;
    ThetaDiff = abs([PitchAngles_buf(1, :); diff(PitchAngles_buf)]);
    for k = 1:N
        buff = cosd(gamma)*Mflap1_store_buf(k)+sind(gamma)*Medge1_store_buf(k);
        max_buff = max(buff);
        max_bufff = max(max_buff, 0);
        PBD = PBD + ThetaDiff(k, 1)*max_bufff^WohlerSlope;
    end
    PitchBearingDamage = PBD / 1e3; % [kNm deg]
    PitchBearningDEL = (PitchBearingDamage ./1e8).^(1/WohlerSlope);
    disp(PitchBearningDEL)
end
function [PBD] = calculatePBD(filter,PitchAngles,Mflap1_store,Medge1_store)
    PitchAngles_buf = PitchAngles(filter:end, :);
    Mflap1_store_buf = Mflap1_store(filter:end, :);
    Medge1_store_buf = Medge1_store(filter:end, :);
    N = length(PitchAngles_buf);
    inverseWohlerSlope = 10;
    PBD = 0;
    gamma = 0:10:360;
    ThetaDiff = [PitchAngles_buf(1, :); diff(PitchAngles_buf)];
    for k = 1:N
        buff = cosd(gamma)*Mflap1_store_buf(k)+sind(gamma)*Medge1_store_buf(k);
        max_buff = max(buff);
        max_bufff = max(max_buff, 0);
        PBD = PBD + ThetaDiff(k, 1)*max_bufff^inverseWohlerSlope;
    end
    PBD = PBD / 1e3;
end
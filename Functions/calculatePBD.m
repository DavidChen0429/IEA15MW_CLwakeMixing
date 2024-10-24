function [PBD] = calculatePBD(PitchAngles,Mflap1_store,Medge1_store)
    N = length(PitchAngles);
    inverseWohlerSlope = 10;
    PBD = 0;
    gamma = 0:10:360;
    ThetaDiff = [PitchAngles(1, :); diff(PitchAngles)];
    for k = 1:N
        buff = cosd(gamma)*Mflap1_store(k)+sind(gamma)*Medge1_store(k);
        max_buff = max(buff);
        max_bufff = max(max_buff, 0);
        PBD = PBD + ThetaDiff(k, 1)*max_bufff^inverseWohlerSlope;
    end
end
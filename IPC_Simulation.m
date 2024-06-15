simFile_loc = [MatlabPath Sim_name_folder '\'];

calllib('QBladeDLL','setLibraryPath',DllPath)
calllib('QBladeDLL','createInstance',2,64)  % NVIDIA GPU
simName = 'NREL5MW_DIPC.sim';
%simName = 'IEA15MW_DIPC.sim';
calllib('QBladeDLL','loadSimDefinition',[simFile_loc simName])
calllib('QBladeDLL','initializeSimulation')

simTime = 1200; % timeStep 0.05
valuestr = 'Pitch Angle Blade 1 [deg]';
valuestr2 = 'Pitch Angle Blade 2 [deg]';
valuestr3 = 'Pitch Angle Blade 3 [deg]';

f = waitbar(0,'Initializing Simulation') ;

for i_for2 = 1:1:simTime
    if mod(i_for2, 100) == 0
        fprintf('%d seconds .\n', i_for2*0.05);
    end

    % Simulation initialization
    calllib('QBladeDLL','advanceTurbineSimulation')

    % Advancess the controller .dll of the selected turbine (num) and 
    % controller outputs are applied to the turbine
    calllib('QBladeDLL','advanceController_at_num',[0 0 0 0 0],0)
    
    % Get the current pitch angle
    Pitch1 = calllib('QBladeDLL','getCustomData_at_num',valuestr, 0, 0) ;
    Pitch2 = calllib('QBladeDLL','getCustomData_at_num',valuestr2, 0, 0); 
    Pitch3 = calllib('QBladeDLL','getCustomData_at_num',valuestr3, 0, 0); 
    
%     % LiDAR data sampling
%     windspeed = ZXTM_lidar(50, 90, 100); % Ring
%      !!! This LIDAR is not the latest version
%     
%     % plot the snapshot
%     x = windspeed.x;
%     y = windspeed.y;
%     z = windspeed.z;
%     speed = windspeed.u_x; % Assuming the field is named 'speed'
%     scatter(y, z, 10, speed, 'filled');
%     xlabel('Y');
%     ylabel('Z');
%     title('LiDAR Wind Speed Visualization');
%     colorbar;
% 
%     % Save the figure to a file
%     filename = sprintf('Data/Figures/Version_ring/figure_%d.png', i_for2);
%     saveas(gcf, filename);

    % Store values in array
    PitchAngles(i_for2,:) = [Pitch1 Pitch2 Pitch3];
    waitbar(i_for2/simTime,f,'Simulation Running')

end
close(f)
calllib('QBladeDLL','closeInstance')
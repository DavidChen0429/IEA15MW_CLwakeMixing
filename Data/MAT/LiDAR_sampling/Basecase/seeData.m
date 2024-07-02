a = load('parallel.mat');
b = load('IEA15_Helix_CCW_Str0.3_U8_Uni_300s_1Dd_1Hz_Circle276_windspeedData.mat');

a = a.LiDAR_data;
b = b.LiDAR_data;

a_ux = arrayfun(@(x) x.u_x(1), a);
b_ux = arrayfun(@(x) x.u_x(1), b);

figure();
plot(a_ux)
hold on
plot(b_ux)

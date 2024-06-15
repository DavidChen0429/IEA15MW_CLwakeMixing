% Define vectors a and b
a = [1; 2; 3]; % Example vector a
b = [4; 5; 6]; % Example vector b

magnitude = norm(b);
b_unit = b/magnitude;

%dot_ab = dot(a, b);
%dot_bb = dot(b, b);
%proj_a_onto_b = (dot_ab / dot_bb) * b;

dot_ab = dot(a, b_unit);
dot_bb = dot(b_unit, b_unit);
proj_a_onto_b = (dot_ab / dot_bb) * b_unit;

figure;

% Plot the original vectors a and b
quiver3(0, 0, 0, a(1), a(2), a(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
hold on;
quiver3(0, 0, 0, b(1), b(2), b(3), 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
quiver3(0, 0, 0, b_unit(1), b_unit(2), b_unit(3), 'LineWidth', 2, 'MaxHeadSize', 0.5);

% Plot the projection of a onto b
quiver3(0, 0, 0, proj_a_onto_b(1), proj_a_onto_b(2), proj_a_onto_b(3), 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);

% Plot the vector from a to its projection on b
line([a(1) proj_a_onto_b(1)], [a(2) proj_a_onto_b(2)], [a(3) proj_a_onto_b(3)], 'Color', 'k', 'LineStyle', '--');

xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Vector Projection Visualization');
legend({'Vector a', 'Vector b','Vector a unit', 'Projection of a onto b', 'Line to Projection'}, 'Location', 'Best');

axis equal;
grid on;
view(3);
hold off;

%% Debug for the LiDAR
x = -100;
num_samples = 200;
y = 0;
z = 105;
x_vel = -8;
y_vel = -2;
z_vel = -1;

pos_original = [x; y; z];
pos = pos_original/norm(pos_original) * 10;
vel = [x_vel; y_vel; z_vel];

% Projection to get u_los
dot_vel_pos = dot(vel, pos);
dot_pos_pos= dot(pos, pos);
u_LOS = (dot_vel_pos / dot_pos_pos) * pos;

% Plot the original vectors a and b
quiver3(0, 0, 0, vel(1), vel(2), vel(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
hold on;
quiver3(0, 0, 0, pos(1), pos(2), pos(3), 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);

% Plot the projection of a onto b
quiver3(0, 0, 0, u_LOS(1), u_LOS(2), u_LOS(3), 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);

% Plot the vector from a to its projection on b
line([vel(1) u_LOS(1)], [vel(2) u_LOS(2)], [vel(3) u_LOS(3)], 'Color', 'k', 'LineStyle', '--');

xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Vector Projection Visualization');
legend({'u_{original}', 'position', 'u_{LOS}', 'Line to Projection'}, 'Location', 'Best');

axis equal;
grid on;
view(3);
hold off;

%% For a 2D scenario
a = [1, 1];
dir = [1, 1];
dir_unit = dir/norm(dir);
(dot(a,dir_unit)/dot(dir_unit,dir_unit))*dir_unit
dot(a, dir_unit)
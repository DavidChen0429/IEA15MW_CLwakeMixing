% This file contains controller design for the decoupled MIMO
% system based on the identified decoupled model sys. This work
% based on the observation that y plays a dominate role in overall ouput

clear
close all
addpath('.\Functions');

% Y dB = 20log10(X) 
% s_c = (2 / Ts) * (z_d - 1) / (z_d + 1);

% Load model
buf_sys = load('Model\RightTransform_Azimuth96\ModelOrder4_noise1p_opposite_decoupled.mat');
A = buf_sys.OLi.A;
B = buf_sys.OLi.B;
C = buf_sys.OLi.C;
D = buf_sys.OLi.D;

sys = buf_sys.OLi;
sys.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
sys.OutputName = {'z_e','y_e'};
G = tf(buf_sys.OLi);        % transfer matrix 
timeStep = 0.1;
DesignOption = 'Matrix2';

% Check result option
showBasicOption = 'N';      % Basic property
showController = 'Y';
showWeight = 'N';
showFreqOption = 'N';       % BD of S, T, U, L
showFreq2Option = 'Y';      % BD of WpS, WuU
showSingularValue = 'N';    % Singular value of S, T
showTimeOption = 'Y';       % Step response
showIteraitve = 'N';
showNyquist = 'N';          % Nyquist stability check
showPerformance = 'N';

% Basic system property
if strcmp(showBasicOption, 'Y')
    fprintf('======== System property \n');
    fprintf(' System Dimension: %.0f \n', size(A, 1));
    fprintf(' Ctrb Matrix Rank: %.0f \n', rank(ctrb(A, B)));
    fprintf(' Obsv Matrix Rank: %.0f \n', rank(obsv(A, C)));
    fprintf(' Eigenvalues of A: \n');
    disp(eig(A));
end

% === RHP zero limitations
% map the discrete zero back to the continuous zero 
zd = 1.0074 + 0.0158i;
sd = (2/timeStep)*(zd-1)/(zd+1);
freq_rad = imag(sd); % Frequency in radians/sec
freq_hz = freq_rad / (2 * pi); % Frequency in Hz

% === Weight function design
% ====== Working weight function
    % W_p = tf([1, 1.6], [100, 1]);  % Emphasizes performance and disturbance rejection
    % W_t = tf([0.01, 1], [6, 1]);   % Emphasizes robustness and noise rejection
    % W_u = tf([1, 1], [50, 1]);     % Emphasizes input magnitude 
% ====== Working weight function
    % Wp0 = tf([1, 0.5], [50, 1]);  % Emphasizes performance and disturbance rejection
    % Wu0 = tf([1, 1], [50, 1]);     % Emphasizes input magnitude 

% 1. Wp(z)
Mp = 10;                % Bound on high freq
Ap = 1.55;               % Bound on low freq
omega_cl = 0.02;        % Closed-loop bandwidth
Wp0 = tf([1/Mp, omega_cl], [1, omega_cl*Ap]);   % Emphasizes performance and disturbance rejection
% 1. Wu(z)
omega_c = 0.5;
B = 10;
Wu0 = 0.4*B^2*tf([1, sqrt(2)*omega_c, omega_c^2], [1, sqrt(2)*B*omega_c, (B*omega_c)^2]);     % Emphasizes input magnitude 
% 3. Wt(z)
Wt0 = tf([1, 1], [5, 1]);   % Emphasizes robustness and noise rejection

close all
% =========== H infinity Control Design
if strcmp(DesignOption, 'Matrix1')  % tried
    Wp_d = c2d(Wp0, timeStep, 'tustin');
    Wt_d = c2d(Wt0, timeStep, 'tustin');
    Wu_d = c2d(Wu0, timeStep, 'tustin');

    P = augw(sys, Wp_d, Wu_d, Wt_d);  % augw creates the weighted augmented plant
    ncont = 2; 
    nmeas = 2; 
    [K,CL,gamma] = hinfsyn(P,nmeas,ncont);
    sys_cl = feedback(sys, K);
elseif strcmp(DesignOption, 'Matrix2')  % less robustness
    Wp_d = c2d(Wp0, timeStep, 'tustin');
    Wu_d = c2d(Wu0, timeStep, 'tustin');
    Wp = blkdiag(Wp_d, Wp_d);
    Wu = blkdiag(Wu_d, Wu_d);
    
    P = augw(sys, Wp, Wu, 0);  % W1:Wp, W2:Wu, W3:Wt.
    ncont = 2; 
    nmeas = 2; 
    [K,CL,gamma] = hinfsyn(P,nmeas,ncont);
    sys_cl = feedback(sys, K);
elseif strcmp(DesignOption, 'Matrix3') % real implementation
    Wp0 = tf([1, 1.6], [100, 1]);
    Wt0 = tf([1, 1], [5, 1]);
    Wp_d = c2d(Wp0, timeStep, 'tustin');
    Wt_d = c2d(Wt0, timeStep, 'tustin');
    Wp = blkdiag(Wp_d, Wp_d);
    Wt = blkdiag(Wt_d, Wt_d);
    
    P = augw(sys, Wp, 0, Wt);  % W1:Wp, W2:Wu, W3:Wt.
    ncont = 2; 
    nmeas = 2; 
    [K,CL,gamma] = hinfsyn(P,nmeas,ncont);
    sys_cl = feedback(sys, K);
end

% =========== Check Result
% Frequency Response
bwG11 = calculateBandwidth(G(1, 1));
bwG22 = calculateBandwidth(G(2, 2));
bw = (bwG11 + bwG22)/2; 
L_mimo = G * K;                     % open-loop tf
S_mimo = feedback(eye(2), L_mimo);  % disturbance --- output
T_mimo = feedback(L_mimo, eye(2));  % reference --- output
U_mimo = feedback(K, G);            % disturbance --- negative control input
L_mimo.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
L_mimo.OutputName = {'z_e','y_e'};
S_mimo.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
S_mimo.OutputName = {'z_e','y_e'};

% Show controller property
if strcmp(showController, 'Y')
    K.InputName = {'e_z', 'e_y'};
    K.OutputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
    figure('Name', 'Controller BD', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bode(K)
    grid on
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)

    figure('Name', 'Controller & U BD', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bodemag(K, U_mimo)
    grid on
    legend('K','U','Location','southeast')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)
end

if strcmp(showWeight, 'Y')
    figure('Name', 'Weight Functions', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bodemag(Wp, Wu);
end

% Frequency Response
if strcmp(showFreqOption, 'Y')
    figure('Name', 'Frequence Result', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bodemag(S_mimo, T_mimo, U_mimo)
    hold on
    axesHandles = findall(gcf, 'Type', 'axes');
    hold off
    grid on
    title('Frequency Response of MIMO System')
    legend('S','T','U','Location','southeast')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)
end

% Check bounded
if strcmp(showFreq2Option, 'Y')
    buf = tf([1/Mp, omega_cl], [1, omega_cl*0.625]);
    buf2 = blkdiag(buf, buf);
    figure('Name', 'Weight Compliance', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bodemag(S_mimo,'m');
    hold on
    bodemag(1/buf2,'m--');
    bodemag(U_mimo,'c');
    bodemag(1/Wu,'c--');
    hold off
    grid on
    legend('S','1/Wp','U','1/Wu','Location','southeast')
    Wp_S_norm = norm(Wp*S_mimo, 'inf');
    Wu_U_norm = norm(Wu*U_mimo, 'inf');
    disp('Nominal Performance:');
    disp(['||Wp * S||_inf = ', num2str(Wp_S_norm)]);
    disp(['||Wu * U||_inf = ', num2str(Wu_U_norm)]);
    
    if Wp_S_norm < 1 && Wu_U_norm < 1
        disp('Nominal Performance: Satisfied.');
    else
        disp('Nominal Performance: NOT satisfied.');
    end
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)
end

% Singular Value
if strcmp(showSingularValue, 'Y')
    figure('Name', 'Singular Value', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    subplot(2, 1, 1)
    sigma(S_mimo);
    hold on
    yline(6,'k--')
    hold off
    ylabel('Magnitude [dB]')
    ylim([-20 10])
    legend('S','6dB','Location','southeast')
    grid on
    subplot(2, 1, 2)
    sigma(T_mimo);
    hold on
    yline(2,'k--')
    hold off
    ylabel('Magnitude [dB]')
    ylim([-100 10])
    title('')
    legend('T','2dB','Location','southeast')
    grid on
    setfigpaper('Width',[30,0.4],'Interpreter','tex','FontSize',20,'linewidth',2)

    % Study the direction 
    [S, wout] = sigma(G); % singular value, rad/s
    condition_number = max(S) ./ min(S);
    bw_rad = bw * 2*pi;
    [~, idx] = min(abs(wout - bw_rad));
    S(:, idx);
    condition_number(idx);
end

% Time Response (Step Response)
if strcmp(showTimeOption, 'Y')
    t = 0:0.1:200;
    figure('Name', 'Step Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    step(sys, t)
    hold on
    step(sys_cl, t)
    axesHandles = findall(gcf, 'Type', 'axes');
    yline(axesHandles(5), 1, 'k--', 'LineWidth', 2);
    yline(axesHandles(4), 0, 'k--', 'LineWidth', 2);
    yline(axesHandles(3), 0, 'k--', 'LineWidth', 2);
    yline(axesHandles(2), 1, 'k--', 'LineWidth', 2);
    hold off
    legend('OL','CL','Location','southeast')
    xlabel('Time [s]')
    ylabel('Magnitude [m]')
%     setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)
end

% Check Stability (Nyquist)
if strcmp(showNyquist, 'Y')
    figure('Name', 'Nyquist Plot', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    Dnyq = eye(size(L_mimo)) + L_mimo;
    nyquist(det(Dnyq));
end

% Check Iterative Time Response
if strcmp(showIteraitve, 'Y')
    A_cl = sys_cl.A;
    B_cl = sys_cl.B;
    C_cl = sys_cl.C;
    D_cl = sys_cl.D;
    t = 0:timeStep:200;
    x = zeros(length(A_cl), length(t)+1);
    y = zeros(length(D_cl), length(t)); 
    r = [5 * ones(1, length(t)+1);
         2 * ones(1, length(t)+1)];
    u = r;
    for i = 1:length(t)
        x(:, i+1) = A_cl * x(:, i) + B_cl * u(:, i);
        y(:, i) = C_cl * x(:, i) + D_cl * u(:, i); 
    end
    figure('Name', 'Iterative Debug', 'NumberTitle', 'off', 'Position', [500, 450, 400, 300]);
    plot(t, y(1,:));  % First output
    hold on;
    plot(t, y(2,:));  % Second output
    hold off
    xlabel('Time (s)');
    ylabel('System Output');
    title('Controlled CL System');
    legend('z_e', 'y_e');
    grid on;
end

%% Analyzing the non-minimum phase zero 
zeros_ij = tzero(G);
nmp_zeros = zeros_ij(abs(zeros_ij) > 1);
for k = 1:length(nmp_zeros)
    zero = nmp_zeros(k);
    omega = angle(zero);
    freq_hz = (omega * 1/timeStep) / (2 * pi);
end

%% Step Response Implementation: Default vs. Iterative
% % close all
% % 1. Default
% closed_loop_sys = feedback(sys, K);
% t = 0:timeStep:200;  % Time vector for simulation
% figure('Name', 'Default Way', 'NumberTitle', 'off', 'Position', [100, 450, 400, 300]);
% [y, tOut] = step(closed_loop_sys, t);
% title('Controlled CL System');
% plot(tOut, y(:, 1) + y(:, 3))
% hold on
% plot(tOut, y(:, 2) + y(:, 4))
% hold off
% xlabel('Time (s)');
% ylabel('System Output');
% title('Controlled CL System');
% legend('z_e', 'y_e');
% grid on;
% 
% % Iterative debug (Below code is right)
% A_cl = closed_loop_sys.A;
% B_cl = closed_loop_sys.B;
% C_cl = closed_loop_sys.C;
% D_cl = closed_loop_sys.D;
% t = 0:timeStep:200;
% x = zeros(length(A_cl), length(t)+1);
% y = zeros(length(D_cl), length(t)); 
% r = [5 * ones(1, length(t)+1);
%      2 * ones(1, length(t)+1)];
% u = r;
% for i = 1:length(t)
%     x(:, i+1) = A_cl * x(:, i) + B_cl * u(:, i);
%     y(:, i) = C_cl * x(:, i) + D_cl * u(:, i); 
% end
% figure('Name', 'Iterative Debug', 'NumberTitle', 'off', 'Position', [500, 450, 400, 300]);
% plot(t, y(1,:));  % First output
% hold on;
% plot(t, y(2,:));  % Second output
% hold off
% xlabel('Time (s)');
% ylabel('System Output');
% title('Controlled CL System');
% legend('z_e', 'y_e');
% grid on;
% 
% % 2. Iterative
% A_K = K.A;
% B_K = K.B;
% C_K = K.C;
% D_K = K.D;
% t = 0:timeStep:200;
% x = zeros(length(A), length(t)+1);
% u = zeros(length(B(1, :)), length(t));
% e = zeros(length(C(:, 1)), length(t));
% y = zeros(length(C(:, 1)), length(t));
% xk = zeros(length(A_K), length(t)+1);
% uk = y; % property of H inf
% yk = zeros(length(C_K(:, 1)), length(t)); % property of H inf  
% r = [5; 2]; 
% for i = 2:length(t)
%     e(:, i) = r - yk(:, i-1);
%     % Update system
%     x(:, i+1) = A * x(:, i) + B * e(:, i);
%     y(:, i) = C * x(:, i) + D * e(:, i); 
%     % Update controller
%     xk(:, i+1) = A_K * xk(:, i) + B_K * y(:, i);
%     yk(:, i) = C_K * xk(:, i) + D_K * y(:, i);
% end
% figure('Name', 'Iterative 1.0', 'NumberTitle', 'off', 'Position', [100, 50, 400, 300]);
% plot(t, y(1,:));  % First output
% hold on;
% plot(t, y(2,:));  % Second output
% plot(t, r(1, :)*ones(1, length(t)), 'k--');
% plot(t, r(2, :)*ones(1, length(t)), 'k--');
% hold off
% xlabel('Time (s)');
% ylabel('System Output');
% title('Controlled CL System');
% legend('z_e', 'y_e');
% grid on;
% 
% % 3. Iterative Alternative
% A_K = K.A;
% B_K = K.B;
% C_K = K.C;
% D_K = K.D;
% t = 0:timeStep:200;
% x = zeros(length(A), length(t)+1);
% u = zeros(length(B(1, :)), length(t));
% e = zeros(length(C(:, 1)), length(t));
% y = zeros(length(C(:, 1)), length(t));
% xk = zeros(length(A_K), length(t)+1);
% uk = y; % property of H inf
% yk = zeros(length(C_K(:, 1)), length(t)); % property of H inf  
% r = [5; 2]; 
% for i = 2:length(t)
%     % Update controller
%     xk(:, i+1) = A_K * xk(:, i) + B_K * y(:, i-1);
%     yk(:, i) = C_K * xk(:, i) + D_K * y(:, i-1);
%     % Get error 
%     e(:, i) = r - yk(:, i);
%     % Update system
%     x(:, i+1) = A * x(:, i) + B * e(:, i);
%     y(:, i) = C * x(:, i) + D * e(:, i); 
% end
% figure('Name', 'Iterative 2.0', 'NumberTitle', 'off', 'Position', [500, 50, 400, 300]);
% plot(t, y(1,:));  % First output
% hold on;
% plot(t, y(2,:));  % Second output
% plot(t, r(1, :)*ones(1, length(t)), 'k--');
% plot(t, r(2, :)*ones(1, length(t)), 'k--');
% hold off
% xlabel('Time (s)');
% ylabel('System Output');
% title('Controlled CL System');
% legend('z_e', 'y_e');
% grid on;
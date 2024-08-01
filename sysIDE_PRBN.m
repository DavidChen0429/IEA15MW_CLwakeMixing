% Parameters
n = 50; 
seedA = [1 0 0 1]; 
seedB = [1 1 0 1]; 
tap = [1 0 0 1]; % polynomial x^4 + x + 1

% PRBN Signal A
stateA = seedA;
prbnA = zeros(1, n);
for i = 1:n
    prbnA(i) = stateA(end); % Output the last bit
    feedbackA = mod(sum(stateA .* tap), 2); % Calculate feedback using tap positions
    stateA = [feedbackA stateA(1:end-1)]; % Shift and insert feedback
end

% PRBN Signal B
stateB = seedB;
prbnB = zeros(1, n);
for i = 1:n
    prbnB(i) = stateB(end); % Output the last bit
    feedbackB = mod(sum(stateB .* tap), 2); % Calculate feedback using tap positions
    stateB = [feedbackB stateB(1:end-1)]; % Shift and insert feedback
end

% Plot the PRBN sequences
figure;
subplot(2, 1, 1);
stairs(prbnA, 'LineWidth', 1.5);
title('Pseudo-Random Binary Noise (PRBN) Sequence A');
xlabel('Sample Index');
ylabel('Binary Value');
ylim([-0.5, 1.5]);
grid on;

subplot(2, 1, 2);
stairs(prbnB, 'LineWidth', 1.5);
title('Pseudo-Random Binary Noise (PRBN) Sequence B');
xlabel('Sample Index');
ylabel('Binary Value');
ylim([-0.5, 1.5]);
grid on;

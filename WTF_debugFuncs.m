timeStep = 0.1;
filter = 3000;

% Debug DEL
neq = 10;
inverseWohlerSlope = 1;
fs = 1/timeStep;
Mflap1_store = Baseline.Mflap1_store;

t = 0:timeStep:10;
x = 2 * sin(pi*t);
plot(x)
cycles_f1 = rainflow(x, fs);
DEL_f1 = (cycles_f1(:,2)'.^inverseWohlerSlope*cycles_f1(:,1)/neq)^(1/inverseWohlerSlope);
disp(DEL_f1)
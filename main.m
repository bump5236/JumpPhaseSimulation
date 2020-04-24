clear
close all
%% setting
set(0, 'DefaultAxesFontSize', 16);
set(0, 'DefaultAxesFontName', 'times');
set(0, 'DefaultTextFontSize', 16);
set(0, 'DefaultTextFontName', 'times');

% add dependent folders
addpath(pwd, 'class');

%% simulation
t0 = 0.0;
t1 = 0.03;
v = 4.2; % [km/h]
h = 0.02;
l1 = 0.07;
l2 = 0.07;
l3 = 0.07;
m1 = 0.2;
m2 = 0.3;
m3 = 0.3;
phi = pi/12;

model = JumpPhase(t0, t1, v, h, l1, l2, l3, m1, m2, m3, phi);

syms t
q1_start = pi/12;
q2_start = pi/6;
q3_start = pi/3;

model.torque(t, q1_start, q2_start, q3_start);

%% plot
fplot(model.tau2, [model.t0 model.t1], 'k', 'LineWidth', 1.2)
xlabel('time [s]', 'Interpreter', 'latex')
ylabel('$\tau$ [Nm]', 'Interpreter', 'latex')
hold on
fplot(model.tau3, [model.t0 model.t1], '--k', 'LineWidth', 1.2)
legend('\tau_{2}', '\tau_{3}', 'Location', 'northwest')
hold off

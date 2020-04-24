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
model = JumpPhase;

syms t
model.vxCoG(t);
model.vzCoG(t);

q1_start = pi/6;
q2_start = pi/6;
q3_start = pi/4;
model.angvel(q1_start, q2_start, q3_start);

model.torque(q1_start, q2_start, q3_start);

%% plot
% tiledlayout(2,1)
% nexttile
% fplot(model.vx, [model.t0 model.t1])
% xlabel('time [s]')
% ylabel('vx [m/s]')
% nexttile
% fplot(model.vz, [model.t0 model.t1])
% xlabel('time [s]')
% ylabel('vz [m/s]')
% nexttile
% fplot(model.dq2, [model.t0 model.t1])
% xlabel('time [s]')
% ylabel('dq2 [rad/s]')
% nexttile
% fplot(model.dq3, [model.t0 model.t1])
% xlabel('time [s]')
% ylabel('dq3 [rad/s]')
% nexttile
% fplot(model.q2, [model.t0 model.t1])
% xlabel('time [s]')
% ylabel('q2 [rad]')
% nexttile
% fplot(model.q3, [model.t0 model.t1])
% xlabel('time [s]')
% ylabel('q3 [rad]')

tiledlayout(3,1)
nexttile
fplot(model.tau1, [model.t0 model.t1])
xlabel('time [s]')
ylabel('\tau1 [Nm]')
nexttile
fplot(model.tau2, [model.t0 model.t1])
xlabel('time [s]')
ylabel('\tau2 [Nm]')
nexttile
fplot(model.tau3, [model.t0 model.t1])
xlabel('time [s]')
ylabel('\tau3 [Nm]')
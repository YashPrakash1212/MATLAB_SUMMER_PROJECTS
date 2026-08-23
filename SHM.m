% Week 5: Computational Physics & Numerical Methods
%  Pendulum Motion involving factors such as Small Angle, Large Angle, Period, Phase Space, Chaos
 
clear; clc; close all;
 
% Important parameters
L  = 1;      % pendulum length (m)
g  = 9.8;    % gravity (m/s^2)
T0 = 2*pi*sqrt(L/g);  % "small angle" period (constant, in theory)
 

% Exercise 1: Simple pendulum, small angle (theta0 = 10 deg)

theta0_1 = 10;      % degrees
dt1      = 0.001;   % s
Tmax1    = 4*T0;    % simulate ~4 periods
 
[theta1, omega1, t1] = simulate_pendulum(theta0_1, 0, L, g, dt1, Tmax1);
 
% Analytical small-angle solution: theta(t) = theta0*cos(sqrt(g/L)*t)
theta0_1_rad = deg2rad(theta0_1);
theta_analytic1 = theta0_1_rad*cos(sqrt(g/L)*t1);
 
figure('Name','Exercise 1: Small-Angle Pendulum');
plot(t1, rad2deg(theta1), 'b-', 'LineWidth', 1.5); hold on;
plot(t1, rad2deg(theta_analytic1), 'r--', 'LineWidth', 1.5); hold off;
grid on;
xlabel('Time (s)'); ylabel('\theta (degrees)');
title(sprintf('Small-Angle Pendulum: \\theta_0 = %d^\\circ', theta0_1));
legend('Euler simulation', 'Analytical (small-angle)', 'Location', 'best');
 
fprintf('--- Exercise 1: Small Angle (%d deg) ---\n', theta0_1);
fprintf('Small-angle period T0 = 2*pi*sqrt(L/g) = %.4f s\n', T0);
fprintf('Simulation vs analytical solution should overlay closely.\n\n');
 
% Exercise 2: Large angle pendulum (theta0 = 170 deg)

theta0_2 = 170;     % since it is 170 degrees the object is basically upside down. 
dt2      = 0.001;
Tmax2    = 20;      % long enough to see the (much longer) real period
 
[theta2, omega2, t2] = simulate_pendulum(theta0_2, 0, L, g, dt2, Tmax2);
 
theta0_2_rad = deg2rad(theta0_2);
theta_analytic2 = theta0_2_rad*cos(sqrt(g/L)*t2);   % small-angle formula (WRONG here on purpose)
 
figure('Name','Exercise 2: Large-Angle Pendulum');
plot(t2, rad2deg(theta2), 'b-', 'LineWidth', 1.5); hold on;
plot(t2, rad2deg(theta_analytic2), 'r--', 'LineWidth', 1.5); hold off;
grid on;
xlabel('Time (s)'); ylabel('\theta (degrees)');
title(sprintf('Large-Angle Pendulum: \\theta_0 = %d^\\circ', theta0_2));
legend('True (Euler) motion', 'Small-angle formula (wrong)', 'Location', 'best');
ylim([-190 190]);
 
Tbig = measure_period(t2, omega2);
 
fprintf('--- Exercise 2: Large Angle (%d deg) ---\n', theta0_2);
fprintf('Small-angle formula predicts T0 = %.4f s (constant, WRONG here)\n', T0);
fprintf('True simulated period       T  = %.4f s\n', Tbig);
fprintf('The real period is noticeably LONGER than T0 -- AP Physics misses this.\n\n');
 
% Exercise 3: Period vs amplitude

theta0_list = [5, 10, 30, 60, 90, 120, 150];  % degrees
nAmp = numel(theta0_list);
dt3  = 0.001;
Tmax3 = 20;   %  window so we always capture a full swing
 
periods_sim   = zeros(nAmp,1);
periods_exact = zeros(nAmp,1);  % elliptic-integral "exact" period, for reference
 
for i = 1:nAmp
    [~, w, tt] = simulate_pendulum(theta0_list(i), 0, L, g, dt3, Tmax3);
    periods_sim(i) = measure_period(tt, w);
 
    % Exact large-angle period via complete elliptic integral of the 1st kind
    th0rad = deg2rad(theta0_list(i));
    k2 = sin(th0rad/2)^2;
    [K, ~] = ellipke(k2);
    periods_exact(i) = 4*sqrt(L/g)*K;
end
 
figure('Name','Exercise 3: Period vs Amplitude');
plot(theta0_list, periods_sim, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
plot(theta0_list, periods_exact, 'g.', 'MarkerSize', 18);
yline(T0, 'r--', 'LineWidth', 1.5);
hold off; grid on;
xlabel('Initial Amplitude \theta_0 (degrees)'); ylabel('Period (s)');
title('Pendulum Period vs Amplitude');
legend('Measured (Euler simulation)', 'Exact (elliptic integral)', ...
       'AP Physics prediction (constant T_0)', 'Location', 'best');
 
fprintf('--- Exercise 3: Period vs Amplitude ---\n');
periodTable = table(theta0_list', periods_sim, periods_exact, ...
    'VariableNames', {'Theta0_deg', 'Period_sim_s', 'Period_exact_s'});
disp(periodTable);
fprintf('AP Physics constant-period prediction: T0 = %.4f s\n', T0);
fprintf('Conclusion: NO -- period grows with amplitude; only ~constant for small angles.\n\n');
 

% Exercise 4: Phase portrait (theta vs omega)

dt4   = 0.001;
Tmax4 = 12;
 
% Oscillating cases (bounded energy -> closed loops), released from rest
osc_theta0_list = [30, 90, 150];   % degrees
 
% Rotating cases (enough kinetic energy to go over the top -> open curves)
% Escape condition (from rest at bottom, theta=0): omega0 > 2*sqrt(g/L)
omega_escape = 2*sqrt(g/L);
rot_omega0_list = [1.3*omega_escape, -1.3*omega_escape];
 
figure('Name','Exercise 4: Phase Portrait'); hold on; grid on;
legendEntries = {};
 
for i = 1:numel(osc_theta0_list)
    [th, w, ~] = simulate_pendulum(osc_theta0_list(i), 0, L, g, dt4, Tmax4);
    plot(rad2deg(th), w, 'LineWidth', 1.5);
    legendEntries{end+1} = sprintf('Oscillating, \\theta_0=%d^\\circ', osc_theta0_list(i)); %#ok<SAGROW>
end
 
for i = 1:numel(rot_omega0_list)
    [th, w, ~] = simulate_pendulum(0, rot_omega0_list(i), L, g, dt4, Tmax4);
    plot(rad2deg(th), w, 'LineWidth', 1.5);
    legendEntries{end+1} = sprintf('Rotating, \\omega_0=%.2f rad/s', rot_omega0_list(i)); %#ok<SAGROW>
end
 
hold off;
xlabel('\theta (degrees)'); ylabel('\omega (rad/s)');
title('Phase Portrait: Oscillation (closed loops) vs Rotation (open curves)');
legend(legendEntries, 'Location', 'best');
 
fprintf('--- Exercise 4: Phase Portrait ---\n');
fprintf('Closed loops = bounded oscillation. Open/unbounded curves = full rotation.\n');
fprintf('Escape angular velocity at theta=0 is omega_escape = %.4f rad/s\n\n', omega_escape);
 
% Exercise 5: Double pendulum (CHAOS!)

m1 = 1; m2 = 1;     % masses (kg)
L1 = 1; L2 = 1;     % lengths (m)
dt5   = 0.0005;     % small dt needed for a stiff, more nonlinear system
Tmax5 = 20;
 
theta1_0 = 120;         % degrees
theta2_0_A = 100;           % run A
theta2_0_B = 100 + 0.001;   % run B: differs by just 0.001 degrees!
 
[t1A, w1A, t2A, w2A, tt5, x2A, y2A] = simulate_double_pendulum(theta1_0, theta2_0_A, m1, m2, L1, L2, g, dt5, Tmax5);
[t1B, w1B, t2B, w2B, ~,   x2B, y2B] = simulate_double_pendulum(theta1_0, theta2_0_B, m1, m2, L1, L2, g, dt5, Tmax5);
 
% Path of the bottom bob, run A
figure('Name','Exercise 5: Double Pendulum Path (Run A)');
plot(x2A, y2A, '-', 'LineWidth', 0.5);
axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)');
title(sprintf('Double Pendulum: Path of Bottom Bob (\\theta_2(0) = %.3f^\\circ)', theta2_0_A));
 
% Overlay both runs to show sensitive dependence on initial conditions
figure('Name','Exercise 5: Chaos - Diverging Paths');
plot(x2A, y2A, 'b-', 'LineWidth', 0.5); hold on;
plot(x2B, y2B, 'r-', 'LineWidth', 0.5); hold off;
axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)');
title('Bottom Bob Path: \theta_2(0) differs by only 0.001^\circ');
legend(sprintf('\\theta_2(0) = %.3f^\\circ', theta2_0_A), ...
       sprintf('\\theta_2(0) = %.3f^\\circ', theta2_0_B), 'Location', 'best');
 
% Quantify divergence over time
divergence = sqrt((x2A - x2B).^2 + (y2A - y2B).^2);
figure('Name','Exercise 5: Divergence vs Time');
semilogy(tt5, divergence, 'k-', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)'); ylabel('Separation distance (m), log scale');
title('Exponential-Looking Divergence of Two Nearly-Identical Double Pendulums');
 
fprintf('--- Exercise 5: Double Pendulum ---\n');
fprintf('Two runs starting only 0.001 deg apart in theta2 end up completely\n');
fprintf('different -- this is the hallmark of chaos (sensitive dependence on\n');
fprintf('initial conditions). Final separation of bottom bob: %.3f m\n\n', divergence(end));
 
 
% Local Functions


function [theta, omega, t] = simulate_pendulum(theta0_deg, omega0, L, g, dt, Tmax)
% SIMULATE_PENDULUM  Semi-implicit ("symplectic") Euler simulation of a simple pendulum:
%  omega_new = omega - (g/L)*sin(theta)*dt
%  theta_new = theta + omega_new*dt
% This form conserves energy better than plain forward Euler, which matters for the period measurements and phase portraits below.
 
    theta0 = deg2rad(theta0_deg);
    N = round(Tmax/dt);
 
    theta = zeros(N+1,1);
    omega = zeros(N+1,1);
    t     = zeros(N+1,1);
 
    theta(1) = theta0;
    omega(1) = omega0;
    t(1)     = 0;
 
    for i = 1:N
        omega(i+1) = omega(i) - (g/L)*sin(theta(i))*dt;
        theta(i+1) = theta(i) + omega(i+1)*dt;
        t(i+1)     = t(i) + dt;
    end
end
 
function T = measure_period(t, omega)
% MEASURE_PERIOD  Estimate the oscillation period from an omega(t) trace.
% For a pendulum released from rest at theta0 (a turning point), omega
% crosses zero (negative -> positive) at exactly the half period. Doubling
% that time gives the full period. Linear interpolation is used between
% samples for sub-timestep accuracy.
 
    idx = find(omega(1:end-1) < 0 & omega(2:end) >= 0, 1, 'first');
    if isempty(idx)
        T = NaN;   % no crossing found (e.g. rotating, not oscillating)
        return;
    end
    t1 = t(idx);   t2 = t(idx+1);
    w1 = omega(idx); w2 = omega(idx+1);
    frac = -w1/(w2 - w1);
    tCross = t1 + frac*(t2 - t1);
    T = 2*tCross;
end
 
function [theta1, omega1, theta2, omega2, t, x2, y2] = simulate_double_pendulum(...
    theta1_0_deg, theta2_0_deg, m1, m2, L1, L2, g, dt, Tmax)
% SIMULATE_DOUBLE_PENDULUM  Euler-method simulation of a double pendulum
% (pendulum 2 hangs from the bob of pendulum 1). Uses the standard
% double-pendulum equations of motion, released from rest.
 
    N = round(Tmax/dt);
 
    theta1 = zeros(N+1,1); omega1 = zeros(N+1,1);
    theta2 = zeros(N+1,1); omega2 = zeros(N+1,1);
    t      = zeros(N+1,1);
 
    theta1(1) = deg2rad(theta1_0_deg);
    theta2(1) = deg2rad(theta2_0_deg);
    omega1(1) = 0;
    omega2(1) = 0;
 
    for i = 1:N
        th1 = theta1(i); th2 = theta2(i);
        w1  = omega1(i); w2  = omega2(i);
 
        denom = 2*m1 + m2 - m2*cos(2*th1 - 2*th2);
 
        num1 = -g*(2*m1 + m2)*sin(th1) - m2*g*sin(th1 - 2*th2) ...
               - 2*sin(th1 - th2)*m2*(w2^2*L2 + w1^2*L1*cos(th1 - th2));
        alpha1 = num1 / (L1*denom);
 
        num2 = 2*sin(th1 - th2)*(w1^2*L1*(m1 + m2) + g*(m1 + m2)*cos(th1) ...
               + w2^2*L2*m2*cos(th1 - th2));
        alpha2 = num2 / (L2*denom);
 
        % Semi-implicit Euler update
        omega1(i+1) = w1 + alpha1*dt;
        omega2(i+1) = w2 + alpha2*dt;
        theta1(i+1) = th1 + omega1(i+1)*dt;
        theta2(i+1) = th2 + omega2(i+1)*dt;
        t(i+1)      = t(i) + dt;
    end
 
    % Cartesian position of the bottom (second) bob
    x1 = L1*sin(theta1);
    y1 = -L1*cos(theta1);
    x2 = x1 + L2*sin(theta2);
    y2 = y1 - L2*cos(theta2);
end
 
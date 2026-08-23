% Week 5: Computational Physics & Numerical Methods
%  Waves on a String: Setup, Propagation, Harmonics, Speed, Interference
%
%  Run the whole file, or step through section-by-section (%% blocks)
%  in the MATLAB Editor / Live Editor.

clear; clc; close all;

% Important parameters
N  = 100;              % number of points along the string
L  = 1;                 % string length (m)
x  = linspace(0, L, N); % positions
dx = x(2) - x(1);

T_tension = 1;   % string tension (N)
mu        = 1;   % linear mass density (kg/m)
c = sqrt(T_tension/mu);   % wave speed c = sqrt(T/mu)

Courant = 0.5;            % Courant number C = c*dt/dx (must be <= 1 for stability)
dt = Courant*dx/c;

fprintf('Wave speed c = sqrt(T/mu) = %.3f m/s\n', c);
fprintf('dx = %.5f m, dt = %.5f s, Courant number = %.3f\n\n', dx, dt, c*dt/dx);

% Exercise 1: Set up the string (plucked at center, triangle shape)

A1 = 0.1;   % pluck amplitude (m)
y0_pluck = triangle_pluck(x, L, A1);

figure('Name','Exercise 1: Initial Condition');
plot(x, y0_pluck, 'b-', 'LineWidth', 2); hold on;
plot([x(1) x(end)], [0 0], 'ko', 'MarkerFaceColor', 'k');
hold off; grid on;
xlabel('x (m)'); ylabel('y (m)');
title('Exercise 1: String Plucked at Center (t = 0)');
ylim([-0.02, 0.13]);

fprintf('--- Exercise 1: Initial Setup ---\n');
fprintf('y(0) = %.4f, y(L) = %.4f (fixed ends), peak at center = %.4f m\n\n', ...
    y0_pluck(1), y0_pluck(end), max(y0_pluck));


% Exercise 2: Simulate wave propagation

v0_zero = zeros(1, N);       % released from rest
Tsim2   = 2.0;                % seconds to simulate (a couple of round trips)
Nsteps2 = round(Tsim2/dt);

[Y2, t2] = simulate_wave(y0_pluck, v0_zero, c, dx, dt, Nsteps2);

frameIdx = round(linspace(1, Nsteps2+1, 6));  % 6 evenly spaced snapshots

figure('Name','Exercise 2: Wave Propagation');
for k = 1:6
    subplot(2,3,k);
    plot(x, Y2(:, frameIdx(k)), 'b-', 'LineWidth', 1.5);
    grid on; ylim([-0.13, 0.13]);
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('t = %.3f s', t2(frameIdx(k))));
end
sgtitle('Exercise 2: Plucked String Evolving in Time');

fprintf('--- Exercise 2: Wave Propagation ---\n');
fprintf('Simulated %.2f s (%d steps). See 6-panel figure for evolution.\n\n', Tsim2, Nsteps2);


% Exercise 3: Standing waves (harmonics)

A3 = 0.1;
figure('Name','Exercise 3: Standing Wave Harmonics');
for n = 1:3
    y0_harm = A3*sin(n*pi*x/L);
    Tn = 2*L/(n*c);              % period of the nth standing-wave mode
    Nsteps_n = round(Tn/dt);

    [Yn, tn] = simulate_wave(y0_harm, v0_zero, c, dx, dt, Nsteps_n);

    subplot(1,3,n); hold on; grid on;
    snapFracs = [0, 0.125, 0.25, 0.375, 0.5];  % fractions of one period
    colors = lines(numel(snapFracs));
    for k = 1:numel(snapFracs)
        idx = round(snapFracs(k)*Nsteps_n) + 1;
        idx = min(idx, Nsteps_n+1);
        plot(x, Yn(:, idx), 'LineWidth', 1.3, 'Color', colors(k,:));
    end
    hold off;
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('Harmonic n = %d (T_%d = %.3f s)', n, n, Tn));
    ylim([-0.13, 0.13]);
    if n == 1
        legend(arrayfun(@(f) sprintf('t = %.2f T', f), snapFracs, 'UniformOutput', false), ...
               'Location', 'south', 'FontSize', 7);
    end
end
sgtitle('Exercise 3: Standing Waves -- Nodes Stay Fixed, Amplitude Oscillates');

fprintf('--- Exercise 3: Standing Wave Harmonics ---\n');
fprintf('Fundamental (n=1), 2nd harmonic (n=2), 3rd harmonic (n=3).\n');
fprintf('Each mode vibrates in place -- nodes never move (standing wave).\n\n');

% Exercise 4: Measure wave speed (traveling pulse)

x0_pulse = 0.15;   % pulse center, away from the boundary
sigma    = 0.03;
Ap4      = 0.1;

y0_pulse4 = Ap4*exp(-((x - x0_pulse)/sigma).^2);
y0_pulse4(1) = 0;                       % keep the fixed boundary exact

dydx = gradient(y0_pulse4, dx);
v0_pulse4 = -c*dydx;                    % initial velocity for a RIGHT-moving pulse only

Nsteps4 = round(1.3*(L - x0_pulse)/c/dt);
[Y4, t4] = simulate_wave(y0_pulse4, v0_pulse4, c, dx, dt, Nsteps4);

% Track the position of the peak over time
xpeak = zeros(1, Nsteps4+1);
for n = 1:Nsteps4+1
    [~, idx] = max(Y4(:,n));
    xpeak(n) = x(idx);
end

% Fit a straight line to xpeak(t) while the pulse is clearly traveling
% (skip the very start, before it's fully separated from the source, and
% the very end, once it nears the far wall)
validIdx = xpeak > (x0_pulse + 0.05) & xpeak < 0.92*L;
p = polyfit(t4(validIdx), xpeak(validIdx), 1);
measured_c = p(1);

figure('Name','Exercise 4: Measuring Wave Speed');
subplot(1,2,1);
plot(x, y0_pulse4, 'b-', 'LineWidth', 1.5); hold on;
plot(x, Y4(:, round(Nsteps4*0.6)), 'r-', 'LineWidth', 1.5);
hold off; grid on;
xlabel('x (m)'); ylabel('y (m)');
title('Pulse: Start vs Later Time'); legend('t = 0', 'later t', 'Location', 'best');

subplot(1,2,2);
plot(t4, xpeak, 'b.'); hold on;
plot(t4(validIdx), polyval(p, t4(validIdx)), 'r-', 'LineWidth', 1.5);
hold off; grid on;
xlabel('Time (s)'); ylabel('Peak position (m)');
title('Peak Position vs Time (slope = measured speed)');
sgtitle('Exercise 4: Wave Speed Measurement');

fprintf('--- Exercise 4: Measuring Wave Speed ---\n');
fprintf('Theoretical c = sqrt(T/mu)      = %.4f m/s\n', c);
fprintf('Measured c (slope of x_peak vs t) = %.4f m/s\n', measured_c);
fprintf('Relative error = %.2f%%\n\n', 100*abs(measured_c-c)/c);


% Exercise 5: Reflection and interference (two pulses collide)

x1 = 0.15;  x2 = 0.85;   % symmetric starting positions
Ap5 = 0.1;

pulseA = Ap5*exp(-((x - x1)/sigma).^2);   % will move right
pulseB = Ap5*exp(-((x - x2)/sigma).^2);   % will move left
y0_combo = pulseA + pulseB;

vA = -c*gradient(pulseA, dx);   % right-moving
vB =  c*gradient(pulseB, dx);   % left-moving
v0_combo = vA + vB;

t_meet = (x2 - x1)/(2*c);       % time for the two pulses to meet at the center
Nsteps5 = round(2.0*t_meet/dt);

[Y5, t5] = simulate_wave(y0_combo, v0_combo, c, dx, dt, Nsteps5);

snapTimes = [0, 0.5*t_meet, t_meet, 1.6*t_meet];
figure('Name','Exercise 5: Pulse Collision');
for k = 1:4
    [~, idx] = min(abs(t5 - snapTimes(k)));
    subplot(2,2,k);
    plot(x, Y5(:, idx), 'b-', 'LineWidth', 1.8);
    grid on; ylim([-0.05, 0.23]);
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('t = %.3f s', t5(idx)));
end
sgtitle('Exercise 5: Two Pulses Pass Through Each Other');

[peakVal, peakStep] = max(max(Y5, [], 1));

fprintf('--- Exercise 5: Reflection and Interference ---\n');
fprintf('Pulses launched from x = %.2f m and x = %.2f m, each amplitude %.2f m.\n', x1, x2, Ap5);
fprintf('Predicted meeting time  t_meet = %.4f s\n', t_meet);
fprintf('Max combined amplitude observed = %.4f m (at t = %.4f s)\n', peakVal, t5(peakStep));
fprintf('Since %.4f m is close to the sum of the two pulse heights (%.2f m),\n', peakVal, 2*Ap5);
fprintf('this confirms constructive interference at the moment of overlap.\n');
fprintf('After overlapping, the two pulses separate again, each still with its\n');
fprintf('original shape and amplitude -- waves pass through each other unchanged.\n\n');


% Local Functions


function y0 = triangle_pluck(x, L, A)
% TRIANGLE_PLUCK  Triangular initial shape: zero at both ends, peak A at
% the center (a string plucked in the middle).
    y0 = zeros(size(x));
    for i = 1:numel(x)
        if x(i) <= L/2
            y0(i) = 2*A*x(i)/L;
        else
            y0(i) = 2*A*(L - x(i))/L;
        end
    end
end

function [Y, t] = simulate_wave(y0, v0, c, dx, dt, Nsteps)
% SIMULATE_WAVE  Explicit finite-difference (leapfrog) solver for the 1D
% wave equation d2y/dt2 = c^2 * d2y/dx2, with fixed (Dirichlet) ends.
%   y(i,n+1) = 2y(i,n) - y(i,n-1) + C^2*(y(i+1,n) - 2y(i,n) + y(i-1,n))
% where C = c*dt/dx is the Courant number (must be <= 1 for stability).
% The first step uses a Taylor-series start-up formula that incorporates
% the initial velocity v0.
%
% Returns Y, an N-by-(Nsteps+1) matrix (each column is a snapshot in
% time), and t, the corresponding time vector.

    N = numel(y0);
    C = c*dt/dx;

    Y = zeros(N, Nsteps+1);
    Y(:,1) = y0(:);

    % Start-up step (uses initial velocity v0)
    Y(2:N-1,2) = Y(2:N-1,1) + dt*v0(2:N-1)' + ...
        0.5*C^2*(Y(3:N,1) - 2*Y(2:N-1,1) + Y(1:N-2,1));
    Y(1,2) = 0; Y(N,2) = 0;   % fixed ends

    % Main leapfrog loop
    for n = 2:Nsteps
        Y(2:N-1,n+1) = 2*Y(2:N-1,n) - Y(2:N-1,n-1) + ...
            C^2*(Y(3:N,n) - 2*Y(2:N-1,n) + Y(1:N-2,n));
        Y(1,n+1) = 0; Y(N,n+1) = 0;
    end

    t = (0:Nsteps)*dt;
end

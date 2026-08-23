% Week 5: Computational Physics & Numerical Methods
%  Orbital Mechanics: Circular/Elliptical Orbits, Kepler's Laws, Escape Velocity

clear; clc; close all;

% Important parameters
G  = 6.674e-11;     % gravitational constant (m^3 kg^-1 s^-2)
M  = 1e30;          % star mass (kg)
GM = G*M;
dt = 1000;           % time step (s), ~16.7 minutes

r0 = [1e11, 0];      % common starting position for Exercises 1, 2, 5 (m)

%% Exercise 1: Circular orbit
v0_1 = [0, 3e4];     % m/s

T1 = orbital_period(r0, v0_1, GM);
N1 = round(1.2*T1/dt);   % simulate ~1.2 orbits

[rx1, ry1, ~, ~, t1] = simulate_orbit(r0, v0_1, GM, dt, N1);

figure('Name','Exercise 1: Circular Orbit');
plot(rx1, ry1, 'b-', 'LineWidth', 1.5); hold on;
plot(0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
hold off; axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)');
title('Exercise 1: Circular Orbit');
legend('Planet orbit', 'Star', 'Location', 'best');

fprintf('--- Exercise 1: Circular Orbit ---\n');
fprintf('Estimated period T = %.3e s (%.1f days)\n\n', T1, T1/86400);

%% Exercise 2: Elliptical orbit
v0_2 = [0, 2e4];     % m/s (slower -> elliptical)

T2 = orbital_period(r0, v0_2, GM);
N2 = round(1.2*T2/dt);

[rx2, ry2, vx2, vy2, t2] = simulate_orbit(r0, v0_2, GM, dt, N2);

figure('Name','Exercise 2: Elliptical Orbit');
plot(rx2, ry2, 'b-', 'LineWidth', 1.5); hold on;
plot(0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
hold off; axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)');
title('Exercise 2: Elliptical Orbit (star at one focus)');
legend('Planet orbit', 'Star (focus)', 'Location', 'best');

fprintf('--- Exercise 2: Elliptical Orbit ---\n');
fprintf('Estimated period T = %.3e s (%.1f days)\n', T2, T2/86400);
fprintf('Notice the star sits at one focus of the ellipse, not the center\n');
fprintf('(confirms Kepler''s 1st law).\n\n');


% Exercise 3: Verify Kepler's 2nd Law (equal areas in equal times)
% Reuse the elliptical orbit from Exercise 2, where speed clearly varies
% around the orbit (fast at perihelion, slow at aphelion).

areaRate = 0.5*(rx2.*vy2 - ry2.*vx2);   % areal velocity dA/dt = 0.5*|r x v|

figure('Name','Exercise 3: Kepler''s 2nd Law');
plot(t2, areaRate, 'b-', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)'); ylabel('Areal velocity, dA/dt (m^2/s)');
title('Exercise 3: Area Swept per Unit Time Should Be Constant');
ylim([0.9*min(areaRate), 1.1*max(areaRate)]);

fprintf('--- Exercise 3: Kepler''s 2nd Law ---\n');
fprintf('Mean areal velocity   = %.4e m^2/s\n', mean(areaRate));
fprintf('Std dev of areal vel. = %.4e m^2/s (%.4f%% of mean)\n', ...
    std(areaRate), 100*std(areaRate)/mean(areaRate));
fprintf('Nearly constant -> confirms equal areas in equal times.\n\n');

% Exercise 4: Multiple planets (Kepler's 3rd law)

rA0 = [1e11, 0];   vA0 = [0, sqrt(GM/norm(rA0))];   % circular orbit, planet A
rB0 = [2e11, 0];   vB0 = [0, sqrt(GM/norm(rB0))];   % circular orbit, planet B (farther)

TA = orbital_period(rA0, vA0, GM);
TB = orbital_period(rB0, vB0, GM);

NA = round(1.05*TA/dt);
NB = round(1.05*TB/dt);

[rxA, ryA, ~, ~, ~] = simulate_orbit(rA0, vA0, GM, dt, NA);
[rxB, ryB, ~, ~, ~] = simulate_orbit(rB0, vB0, GM, dt, NB);

figure('Name','Exercise 4: Multiple Planets');
plot(rxA, ryA, 'b-', 'LineWidth', 1.5); hold on;
plot(rxB, ryB, 'r-', 'LineWidth', 1.5);
plot(0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
hold off; axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)');
title('Exercise 4: Two Planets Orbiting the Same Star');
legend('Planet A (r = 1e11 m)', 'Planet B (r = 2e11 m)', 'Star', 'Location', 'best');

keplerConstA = TA^2/norm(rA0)^3;
keplerConstB = TB^2/norm(rB0)^3;
theoreticalConst = 4*pi^2/GM;

fprintf('--- Exercise 4: Kepler''s 3rd Law ---\n');
fprintf('Planet A: r = %.2e m, T = %.3e s (%.1f days)\n', norm(rA0), TA, TA/86400);
fprintf('Planet B: r = %.2e m, T = %.3e s (%.1f days)\n', norm(rB0), TB, TB/86400);
fprintf('Farther planet (B) has the longer period, as expected.\n');
fprintf('T^2/r^3 for A = %.4e s^2/m^3\n', keplerConstA);
fprintf('T^2/r^3 for B = %.4e s^2/m^3\n', keplerConstB);
fprintf('Theoretical 4*pi^2/GM = %.4e s^2/m^3  (all three should match)\n\n', theoreticalConst);


% Exercise 5: Escape velocity

r_start = norm(r0);
v_escape = sqrt(2*GM/r_start);

v_below = 0.95*v_escape;   % should stay bound (returns)
v_above = 1.05*v_escape;   % should escape (never returns)

v0_below = [0, v_below];
v0_above = [0, v_above];

T_below = orbital_period(r0, v0_below, GM);   % finite (bound ellipse)
Tmax5   = 1.1*T_below;                        % same duration used for both, for fair comparison
N5      = round(Tmax5/dt);

[rx_b, ry_b, ~, ~, t5] = simulate_orbit(r0, v0_below, GM, dt, N5);
[rx_a, ry_a, ~, ~, ~ ] = simulate_orbit(r0, v0_above, GM, dt, N5);

dist_below = sqrt(rx_b.^2 + ry_b.^2);
dist_above = sqrt(rx_a.^2 + ry_a.^2);

figure('Name','Exercise 5: Escape Velocity - Trajectories');
plot(rx_b, ry_b, 'b-', 'LineWidth', 1.5); hold on;
plot(rx_a, ry_a, 'r-', 'LineWidth', 1.5);
plot(0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
hold off; axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)');
title('Exercise 5: Below vs Above Escape Velocity');
legend(sprintf('v = %.0f m/s (95%% v_{esc}) -- bound', v_below), ...
       sprintf('v = %.0f m/s (105%% v_{esc}) -- escapes', v_above), ...
       'Star', 'Location', 'best');

figure('Name','Exercise 5: Escape Velocity - Distance vs Time');
plot(t5, dist_below, 'b-', 'LineWidth', 1.5); hold on;
plot(t5, dist_above, 'r-', 'LineWidth', 1.5);
hold off; grid on;
xlabel('Time (s)'); ylabel('Distance from star (m)');
title('Exercise 5: Distance from Star Over Time');
legend('Below escape velocity (returns)', 'Above escape velocity (flies away)', 'Location', 'best');

fprintf('--- Exercise 5: Escape Velocity ---\n');
fprintf('Theoretical escape velocity v_esc = sqrt(2*GM/r) = %.1f m/s\n', v_escape);
fprintf('Below-escape run: v = %.1f m/s -> bound elliptical orbit, distance oscillates (returns).\n', v_below);
fprintf('Above-escape run: v = %.1f m/s -> unbound trajectory, distance grows without bound (escapes).\n\n', v_above);



% Local Functions


function [rx, ry, vx, vy, t] = simulate_orbit(r0, v0, GM, dt, N)
% SIMULATE_ORBIT  Semi-implicit ("symplectic") Euler simulation of a
% planet orbiting a fixed star at the origin.
%   a = -GM*r/|r|^3
%   v_new = v + a*dt
%   r_new = r + v_new*dt
% Using the just-updated velocity to update position (rather than the old
% velocity) keeps the orbit's energy nearly conserved, which is essential
% for realistic-looking closed orbits over many steps.

    rx = zeros(N+1,1); ry = zeros(N+1,1);
    vx = zeros(N+1,1); vy = zeros(N+1,1);
    t  = zeros(N+1,1);

    rx(1) = r0(1); ry(1) = r0(2);
    vx(1) = v0(1); vy(1) = v0(2);
    t(1)  = 0;

    for i = 1:N
        rMag = sqrt(rx(i)^2 + ry(i)^2);
        ax = -GM*rx(i)/rMag^3;
        ay = -GM*ry(i)/rMag^3;

        vx(i+1) = vx(i) + ax*dt;
        vy(i+1) = vy(i) + ay*dt;
        rx(i+1) = rx(i) + vx(i+1)*dt;
        ry(i+1) = ry(i) + vy(i+1)*dt;
        t(i+1)  = t(i) + dt;
    end
end

function T = orbital_period(r0, v0, GM)
% ORBITAL_PERIOD  Period of a two-body orbit via the vis-viva equation
% (to get the semi-major axis) and Kepler's 3rd law. Returns Inf for
% unbound (parabolic/hyperbolic) trajectories.

    r = norm(r0);
    v = norm(v0);
    invA = 2/r - v^2/GM;   % 1/a from vis-viva
    if invA <= 0
        T = Inf;
    else
        a = 1/invA;
        T = 2*pi*sqrt(a^3/GM);
    end
end
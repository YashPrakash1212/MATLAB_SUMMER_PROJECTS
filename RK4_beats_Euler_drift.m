% Task 3 & 4: RK4 vs Euler for Orbital Motion, and Energy Conservation
%  Task 3: Implement RK4 for the orbit simulation and compare against Euler's method : Euler drifts, RK4 stays closed.
%  Task 4: Track KE + PE over time for both methods : Euler's total energy drifts, RK4's stays essentially constant.
%  Run the whole file, or step through section-by-section (%% blocks)
%  in the MATLAB Editor / Live Editor.
 
clear; clc; close all;
 
% Shared parameters (same orbit as before) 
G  = 6.674e-11;     % gravitational constant (m^3 kg^-1 s^-2)
M  = 1e30;          % star mass (kg)
GM = G*M;
m_planet = 1;        % planet mass (kg) -- arbitrary; only affects absolute
                      % energy scale, not the qualitative drift-vs-conserved story
 
r0 = [1e11, 0];       % starting position (m)
v0 = [0, 3e4];        % starting velocity (m/s), near-circular
dt = 1000;             % time step (s)
 
T_orbit = orbital_period(r0, v0, GM);
nOrbits = 5;
N = round(nOrbits*T_orbit/dt);   % simulate several full orbits so drift is obvious
 
fprintf('Orbital period T = %.3e s (%.1f days)\n', T_orbit, T_orbit/86400);
fprintf('Simulating %d orbits: %d steps at dt = %d s\n\n', nOrbits, N, dt);
 
% Task 3: Implement RK4, compare against (naive) Euler

[rxE, ryE, vxE, vyE, tE] = simulate_orbit_euler(r0, v0, GM, dt, N);
[rxR, ryR, vxR, vyR, tR] = simulate_orbit_rk4(r0, v0, GM, dt, N);
 
figure('Name','Task 3: Euler vs RK4 Orbit Shape');
plot(rxE, ryE, 'r-', 'LineWidth', 1); hold on;
plot(rxR, ryR, 'b-', 'LineWidth', 1);
plot(0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
hold off; axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)');
title(sprintf('Task 3: Orbit Shape over %d Orbits -- Euler Drifts, RK4 Stays Closed', nOrbits));
legend('Euler (naive)', 'RK4', 'Star', 'Location', 'best');
 
% Quantify the drift: how far has each method's final radius moved from the start?
r_start = norm(r0);
r_final_euler = norm([rxE(end), ryE(end)]);
r_final_rk4   = norm([rxR(end), ryR(end)]);
 
fprintf('--- Task 3: Euler vs RK4 ---\n');
fprintf('Starting radius            = %.4e m\n', r_start);
fprintf('Euler radius after %d orbits = %.4e m (drifted by %.2f%%)\n', ...
    nOrbits, r_final_euler, 100*(r_final_euler-r_start)/r_start);
fprintf('RK4 radius after %d orbits   = %.4e m (drifted by %.2f%%)\n\n', ...
    nOrbits, r_final_rk4, 100*(r_final_rk4-r_start)/r_start);
 

% Task 4: Energy conservation check

E_euler = compute_energy(rxE, ryE, vxE, vyE, GM, m_planet);
E_rk4   = compute_energy(rxR, ryR, vxR, vyR, GM, m_planet);
 
figure('Name','Task 4: Energy Conservation');
plot(tE/86400, E_euler, 'r-', 'LineWidth', 1.3); hold on;
plot(tR/86400, E_rk4,   'b-', 'LineWidth', 1.3);
hold off; grid on;
xlabel('Time (days)'); ylabel('Total Energy, KE + PE (J)');
title('Task 4: Total Orbital Energy vs Time');
legend('Euler (naive) -- drifts', 'RK4 -- nearly constant', 'Location', 'best');
 
E0_euler = E_euler(1);
E0_rk4   = E_rk4(1);
drift_euler_pct = 100*(E_euler(end) - E0_euler)/abs(E0_euler);
drift_rk4_pct   = 100*(E_rk4(end)   - E0_rk4)/abs(E0_rk4);
 
fprintf('--- Task 4: Energy Conservation ---\n');
fprintf('Initial energy            = %.4e J\n', E0_euler);
fprintf('Euler: final energy       = %.4e J  (drift = %.3f%%)\n', E_euler(end), drift_euler_pct);
fprintf('RK4:   final energy       = %.4e J  (drift = %.5f%%)\n', E_rk4(end), drift_rk4_pct);
fprintf('\nConclusion: Euler''s error accumulates step after step, steadily\n');
fprintf('pumping energy into (or draining it from) the orbit, so the planet\n');
fprintf('spirals in or out. RK4''s much smaller local error means the energy\n');
fprintf('stays essentially flat and the orbit stays closed -- this is exactly\n');
fprintf('why higher-order integrators matter for long-term simulations.\n\n');
 
 

% Local Functions

 
function [rx, ry, vx, vy, t] = simulate_orbit_euler(r0, v0, GM, dt, N)
% SIMULATE_ORBIT_EULER  Naive (fully explicit) Euler's method:
%   a       = -GM*r/|r|^3          (computed from the CURRENT state)
%   r_new   = r + v*dt             (uses the OLD velocity)
%   v_new   = v + a*dt
% This is the textbook Euler's method, and it systematically injects
% energy into the orbit, causing it to spiral outward over time.
 
    rx = zeros(N+1,1); ry = zeros(N+1,1);
    vx = zeros(N+1,1); vy = zeros(N+1,1);
    t  = zeros(N+1,1);
 
    rx(1) = r0(1); ry(1) = r0(2);
    vx(1) = v0(1); vy(1) = v0(2);
 
    for i = 1:N
        rMag = sqrt(rx(i)^2 + ry(i)^2);
        ax = -GM*rx(i)/rMag^3;
        ay = -GM*ry(i)/rMag^3;
 
        rx(i+1) = rx(i) + vx(i)*dt;   % old velocity used here
        ry(i+1) = ry(i) + vy(i)*dt;
        vx(i+1) = vx(i) + ax*dt;
        vy(i+1) = vy(i) + ay*dt;
        t(i+1)  = t(i) + dt;
    end
end
 
function [rx, ry, vx, vy, t] = simulate_orbit_rk4(r0, v0, GM, dt, N)
% SIMULATE_ORBIT_RK4  Classical 4th-order Runge-Kutta integration of the
% orbital equations of motion, state = [x, y, vx, vy].
 
    rx = zeros(N+1,1); ry = zeros(N+1,1);
    vx = zeros(N+1,1); vy = zeros(N+1,1);
    t  = zeros(N+1,1);
 
    state = [r0(1); r0(2); v0(1); v0(2)];
    rx(1) = state(1); ry(1) = state(2);
    vx(1) = state(3); vy(1) = state(4);
 
    for i = 1:N
        k1 = orbit_deriv(state,          GM);
        k2 = orbit_deriv(state + dt/2*k1, GM);
        k3 = orbit_deriv(state + dt/2*k2, GM);
        k4 = orbit_deriv(state + dt*k3,   GM);
 
        state = state + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
 
        rx(i+1) = state(1); ry(i+1) = state(2);
        vx(i+1) = state(3); vy(i+1) = state(4);
        t(i+1)  = t(i) + dt;
    end
end
 
function dstate = orbit_deriv(state, GM)
% ORBIT_DERIV  Right-hand side of the orbital ODE system:
%   dx/dt = vx,  dy/dt = vy,  dvx/dt = ax,  dvy/dt = ay
    x = state(1); y = state(2);
    vx = state(3); vy = state(4);
 
    rMag = sqrt(x^2 + y^2);
    ax = -GM*x/rMag^3;
    ay = -GM*y/rMag^3;
 
    dstate = [vx; vy; ax; ay];
end
 
function E = compute_energy(rx, ry, vx, vy, GM, m)
% COMPUTE_ENERGY  Total mechanical energy (KE + PE) of an orbiting body.
%   KE = 0.5*m*v^2
%   PE = -G*M*m/r = -GM*m/r
    r  = sqrt(rx.^2 + ry.^2);
    v2 = vx.^2 + vy.^2;
    KE = 0.5*m*v2;
    PE = -GM*m./r;
    E  = KE + PE;
end
 
function T = orbital_period(r0, v0, GM)
% ORBITAL_PERIOD  Period of a two-body orbit via the vis-viva equation
% (semi-major axis) and Kepler's 3rd law.
    r = norm(r0);
    v = norm(v0);
    invA = 2/r - v^2/GM;
    if invA <= 0
        T = Inf;
    else
        a = 1/invA;
        T = 2*pi*sqrt(a^3/GM);
    end
end
 
% Week 5: Computational Physics & Numerical Methods
%  Projectile Motion: With and Without Air Resistance

clear; clc; close all;

% Important parameters
g       = 9.81;    % gravity (in m/s^2)
v0      = 50;       % initial speed (in m/s)
m       = 1;        % projectile mass (in kg) -- b is defined per unit mass
dt      = 0.01;      % Euler time step (in seconds)

% Exercise 1: Projectile WITHOUT air resistance
theta1_deg = 45;
theta1     = deg2rad(theta1_deg);

Tflight1 = 2*v0*sin(theta1)/g;              % total flight time
Range1   = v0^2*sin(2*theta1)/g;            % horizontal range
MaxH1    = (v0*sin(theta1))^2/(2*g);        % max height

t1 = linspace(0, Tflight1, 1000);
x1 = v0*cos(theta1)*t1;
y1 = v0*sin(theta1)*t1 - 0.5*g*t1.^2;

figure('Name','Exercise 1: No Air Resistance');
plot(x1, y1, 'b-', 'LineWidth', 2); grid on;
xlabel('Horizontal Distance (m)'); ylabel('Height (m)');
title(sprintf('Projectile Motion - No Air Resistance (\\theta = %d^\\circ)', theta1_deg));
axis equal;

fprintf('--- Exercise 1: No Air Resistance ---\n');
fprintf('Range      = %.2f m\n', Range1);
fprintf('Max Height = %.2f m\n', MaxH1);
fprintf('Flight Time= %.2f s\n\n', Tflight1);

% Exercise 2: Projectile WITH air resistance (Euler's method)
b2 = 0.1;   % drag coefficient

[x2, y2, t2, Range2, MaxH2, Tflight2] = simulate_with_air(v0, theta1_deg, b2, m, g, dt);

figure('Name','Exercise 2: With Air Resistance');
plot(x2, y2, 'r-', 'LineWidth', 2); grid on;
xlabel('Horizontal Distance (m)'); ylabel('Height (m)');
title(sprintf('Projectile Motion - With Air Resistance (b = %.2f)', b2));

fprintf('--- Exercise 2: With Air Resistance (b = %.2f) ---\n', b2);
fprintf('Range      = %.2f m\n', Range2);
fprintf('Max Height = %.2f m\n', MaxH2);
fprintf('Flight Time= %.2f s\n\n', Tflight2);

% Exercise 3: Compare both trajectories on the same plot

figure('Name','Exercise 3: Comparison');
plot(x1, y1, 'b-', 'LineWidth', 2); hold on;
plot(x2, y2, 'r--', 'LineWidth', 2); hold off;
grid on;
xlabel('Horizontal Distance (m)'); ylabel('Height (m)');
title('Projectile Motion: No Air vs With Air');
legend('No Air', 'With Air', 'Location', 'best');

rangeDiff    = Range1 - Range2;
rangeDiffPct = 100 * rangeDiff / Range1;

fprintf('--- Exercise 3: Comparison ---\n');
fprintf('Range without air : %.2f m\n', Range1);
fprintf('Range with air    : %.2f m\n', Range2);
fprintf('Range reduction   : %.2f m (%.1f%% shorter)\n\n', rangeDiff, rangeDiffPct);

%% Exercise 4: Vary the drag coefficient
b_vals = [0, 0.05, 0.1, 0.2, 0.5];
nB     = numel(b_vals);

trajX      = cell(nB,1);
trajY      = cell(nB,1);
RangeVals  = zeros(nB,1);
MaxHVals   = zeros(nB,1);
TimeVals   = zeros(nB,1);

for i = 1:nB
    [xi, yi, ~, Ri, Hi, Ti] = simulate_with_air(v0, theta1_deg, b_vals(i), m, g, dt);
    trajX{i}     = xi;
    trajY{i}     = yi;
    RangeVals(i) = Ri;
    MaxHVals(i)  = Hi;
    TimeVals(i)  = Ti;
end

figure('Name','Exercise 4: Varying Drag Coefficient');
hold on; grid on;
colors = lines(nB);
legendEntries = strings(nB,1);
for i = 1:nB
    plot(trajX{i}, trajY{i}, 'LineWidth', 2, 'Color', colors(i,:));
    legendEntries(i) = sprintf('b = %.2f', b_vals(i));
end
hold off;
xlabel('Horizontal Distance (m)'); ylabel('Height (m)');
title('Projectile Trajectories for Different Drag Coefficients');
legend(legendEntries, 'Location', 'best');

fprintf('--- Exercise 4: Varying Drag Coefficient ---\n');
dragTable = table(b_vals', RangeVals, MaxHVals, TimeVals, ...
    'VariableNames', {'b', 'Range_m', 'MaxHeight_m', 'FlightTime_s'});
disp(dragTable);
fprintf('\n');

% Exercise 5: Find optimal launch angle WITH air resistance

b5     = 0.1;
angles = 20:1:70;
nA     = numel(angles);
Range5 = zeros(nA,1);

for i = 1:nA
    [~, ~, ~, Ri, ~, ~] = simulate_with_air(v0, angles(i), b5, m, g, dt);
    Range5(i) = Ri;
end

[bestRange, idx] = max(Range5);
bestAngle = angles(idx);

figure('Name','Exercise 5: Optimal Launch Angle');
plot(angles, Range5, 'b-', 'LineWidth', 2); hold on;
plot(bestAngle, bestRange, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on;
xlabel('Launch Angle (degrees)'); ylabel('Range (m)');
title(sprintf('Range vs Launch Angle (with air, b = %.2f)', b5));
legend('Range', sprintf('Optimal: %d^\\circ', bestAngle), 'Location', 'best');
hold off;

fprintf('--- Exercise 5: Optimal Launch Angle (With Air) ---\n');
fprintf('Optimal angle = %d degrees\n', bestAngle);
fprintf('Max range     = %.2f m\n', bestRange);
fprintf('(Without air, optimal angle is always 45 degrees.)\n\n');


%Local Functions

function [x, y, t, R, maxH, Tflight] = simulate_with_air(v0, theta_deg, b, m, g, dt)
% Projectile simulation with linear drag, using Euler Method
%   Fx = -b*vx, Fy = -m*g - b*vy
%   [x,y,t,R,maxH,Tflight] = simulate_with_air(v0,theta_deg,b,m,g,dt)

    theta = deg2rad(theta_deg);

    % initial state
    vx = v0*cos(theta);
    vy = v0*sin(theta);
    x  = 0;
    y  = 0;
    t  = 0;

    xs = x; ys = y; ts = t;

    % Euler integration loop: while y >= 0, update ax, ay, vx, vy, x, y
    while y >= 0
        ax = -b*vx/m;
        ay = -g - b*vy/m;

        vx = vx + ax*dt;
        vy = vy + ay*dt;
        x  = x  + vx*dt;
        y  = y  + vy*dt;
        t  = t  + dt;

        xs(end+1) = x; %#ok<AGROW>
        ys(end+1) = y; %#ok<AGROW>
        ts(end+1) = t; %#ok<AGROW>
    end

    % Linear interpolation between the last two points to land exactly at y = 0
    x1 = xs(end-1); y1 = ys(end-1); t1 = ts(end-1);
    x2 = xs(end);   y2 = ys(end);   t2 = ts(end);
    frac = y1 / (y1 - y2);

    R       = x1 + frac*(x2 - x1);
    Tflight = t1 + frac*(t2 - t1);
    maxH    = max(ys);

    % clean up the final point so the plotted trajectory ends at y = 0
    xs(end) = R;
    ys(end) = 0;
    ts(end) = Tflight;

    x = xs; y = ys; t = ts;
end

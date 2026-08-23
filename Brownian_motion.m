% Exercise 1: 2D random walk (Brownian motion)
rng(0);
n_steps = 10000;

angles = 2*pi*rand(n_steps, 1);
step_sizes = rand(n_steps, 1);

dx = step_sizes .* cos(angles);
dy = step_sizes .* sin(angles);

x = [0; cumsum(dx)];
y = [0; cumsum(dy)];

figure;
plot(x, y, 'LineWidth', 0.5, 'Color', [0.2 0.4 0.8]); hold on;
plot(0, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(x(end), y(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('x');
ylabel('y');
title(sprintf('2D Random Walk (%d steps) - Brownian Motion', n_steps));
legend('Path', 'Start', 'End');
axis equal;
grid on;




% Exercise 2: Diffusion of many molecules
n_molecules = 1000;
n_steps = 2000;
checkpoints = [1, 101, 501, 2001];  % t=0,100,500,2000 (1-indexed)

angles = 2*pi*rand(n_steps, n_molecules);
steps = rand(n_steps, n_molecules);
dx = steps .* cos(angles);
dy = steps .* sin(angles);

X = [zeros(1, n_molecules); cumsum(dx, 1)];
Y = [zeros(1, n_molecules); cumsum(dy, 1)];

figure;
titles = {'t = 0', 't = 100', 't = 500', 't = 2000'};
for k = 1:4
    subplot(1, 4, k);
    scatter(X(checkpoints(k), :), Y(checkpoints(k), :), 8, [0.9 0.4 0.1], 'filled');
    title(titles{k});
    xlim([-60 60]); ylim([-60 60]);
    axis square;
    grid on;
end
sgtitle('Diffusion of 1000 Molecules Over Time');






% Exercise 3: Einstein relation -- <x^2> vs t
times = [10, 50, 100, 200, 400, 800, 1600];
n_walks = 1000;
msd = zeros(size(times));

for k = 1:length(times)
    t = times(k);
    angles = 2*pi*rand(t, n_walks);
    steps = rand(t, n_walks);
    dx = steps .* cos(angles);
    dy = steps .* sin(angles);
    x_final = sum(dx, 1);
    y_final = sum(dy, 1);
    sq_disp = x_final.^2 + y_final.^2;
    msd(k) = mean(sq_disp);
end

p = polyfit(times, msd, 1);
slope = p(1);
intercept = p(2);
D = slope / 2;  % slope = 2D in 2D

figure;
plot(times, msd, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
plot(times, slope*times + intercept, 'r--', 'LineWidth', 2);
xlabel('Time (steps)');
ylabel('Mean Squared Displacement <x^2>');
title('Einstein Relation: <x^2> vs Time');
legend('Simulated', sprintf('Fit: slope = %.3f', slope), 'Location', 'northwest');
grid on;

fprintf('Slope = %.4f, Diffusion coefficient D = %.4f\n', slope, D);





% Exercise 4: Diffusion PDE -- spike flattening into a Gaussian
D = 0.5;
L = 100;
dx = 1.0;
dt = 0.2;
nx = round(L/dx);

C = zeros(nx, 1);
C(round(nx/2)) = 1000;  % delta function spike at center

alpha = D*dt/dx^2;
assert(alpha < 0.5, 'Unstable! Reduce dt or increase dx.');

save_steps = [1, 51, 201, 501];
snapshots = zeros(nx, length(save_steps));
save_idx = 1;

for step = 1:501
    if any(step == save_steps)
        snapshots(:, save_idx) = C;
        save_idx = save_idx + 1;
    end
    C_new = C;
    C_new(2:end-1) = C(2:end-1) + alpha*(C(3:end) - 2*C(2:end-1) + C(1:end-2));
    C = C_new;
end

figure; hold on;
labels = {'t=0', 't=50', 't=200', 't=500'};
for k = 1:size(snapshots, 2)
    plot(snapshots(:, k), 'LineWidth', 2, 'DisplayName', labels{k});
end
xlabel('Position');
ylabel('Concentration');
title('Diffusion: Spike Flattening into a Gaussian');
legend show;
grid on;





% Exercise 5: Diffusion across a membrane
k = 0.05;        % membrane permeability
C_left = 100.0;
C_right = 0.0;
n_t = 500;

left_hist = zeros(n_t+1, 1);
right_hist = zeros(n_t+1, 1);
left_hist(1) = C_left;
right_hist(1) = C_right;

for t = 1:n_t
    flux = k * (C_left - C_right);
    C_left = C_left - flux;
    C_right = C_right + flux;
    left_hist(t+1) = C_left;
    right_hist(t+1) = C_right;
end

figure;
plot(left_hist, 'LineWidth', 2); hold on;
plot(right_hist, 'LineWidth', 2);
xlabel('Time');
ylabel('Concentration');
title("Diffusion Across a Membrane (Fick's Law)");
legend('Left compartment (e.g. blood)', 'Right compartment (e.g. alveoli)');
grid on;





% Exercise 6: Stock (Brownian motion + drift) vs molecule (no drift)
n_steps = 1000;

% Stock: Brownian motion WITH drift
drift = 0.001;
volatility = 1.0;
stock_returns = drift + volatility * randn(n_steps, 1);
stock_path = [100; 100 + cumsum(stock_returns)];

% Molecule: Brownian motion WITHOUT drift
molecule_steps = randn(n_steps, 1);
molecule_path = [0; cumsum(molecule_steps)];

figure;
subplot(1, 2, 1);
plot(stock_path, 'g-', 'LineWidth', 1.5);
title('Stock Price (Brownian motion + drift)');
xlabel('Time');
ylabel('Price');
grid on;

subplot(1, 2, 2);
plot(molecule_path, 'b-', 'LineWidth', 1.5);
title('Molecule Position (Brownian motion, no drift)');
xlabel('Time');
ylabel('Position');
grid on;
% Exercise 1: Exponential growth
r = 0.5;
N0 = 100;
dt = 0.01;
t_end = 20;
t = 0:dt:t_end;
N = zeros(size(t));
N(1) = N0;

for i = 1:length(t)-1
    dN = r * N(i);
    N(i+1) = N(i) + dN*dt;
end

figure;
plot(t, N, 'b-', 'LineWidth', 2);
xlabel('Time');
ylabel('Population N');
title('Exponential Growth: dN/dt = rN');
grid on;

% Find when N hits 1 million
idx = find(N >= 1e6, 1);
if ~isempty(idx)
    fprintf('Population hits 1,000,000 at t = %.2f\n', t(idx));
else
    t_analytical = log(1e6/N0) / r;
    fprintf('Analytical: hits 1,000,000 at t = %.2f (not reached in simulation window)\n', t_analytical);
end



% Exercise 2: Logistic growth vs exponential
r = 0.5;
N0 = 100;
K = 10000;
dt = 0.01;
t_end = 30;
t = 0:dt:t_end;

N_exp = zeros(size(t));
N_log = zeros(size(t));
N_exp(1) = N0;
N_log(1) = N0;

for i = 1:length(t)-1
    dN_exp = r * N_exp(i);
    dN_log = r * N_log(i) * (1 - N_log(i)/K);
    N_exp(i+1) = N_exp(i) + dN_exp*dt;
    N_log(i+1) = N_log(i) + dN_log*dt;
end

% Inflection point occurs where N = K/2
[~, idx_inflect] = min(abs(N_log - K/2));
t_inflect = t(idx_inflect);

figure;
plot(t, N_exp, 'r--', 'LineWidth', 2); hold on;
plot(t, N_log, 'b-', 'LineWidth', 2);
plot(t_inflect, K/2, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
yline(K, 'g--', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Population N');
title('Exponential vs Logistic Growth');
legend('Exponential', 'Logistic', 'Inflection point', 'Carrying capacity K', 'Location', 'northwest');
grid on;

fprintf('Inflection point at t = %.2f, N = %.0f (half of K)\n', t_inflect, K/2);


% Exercise 3: Lotka-Volterra predator-prey
a = 1.0; b = 0.1;      % prey growth, predation rate
c = 1.5; d = 0.075;    % predator death, predator growth from eating prey
dt = 0.01;
t_end = 40;
t = 0:dt:t_end;

prey = zeros(size(t));
pred = zeros(size(t));
prey(1) = 40;
pred(1) = 9;

for i = 1:length(t)-1
    dprey = a*prey(i) - b*prey(i)*pred(i);
    dpred = -c*pred(i) + d*prey(i)*pred(i);
    prey(i+1) = prey(i) + dprey*dt;
    pred(i+1) = pred(i) + dpred*dt;
end

figure;
plot(t, prey, 'b-', 'LineWidth', 2); hold on;
plot(t, pred, 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Population');
title('Lotka-Volterra Predator-Prey Dynamics');
legend('Prey (rabbits)', 'Predators (foxes)');
grid on;



% Exercise 4: Phase portrait of predator-prey
a = 1.0; b = 0.1; c = 1.5; d = 0.075;
dt = 0.01;
t_end = 40;
t = 0:dt:t_end;

initial_conditions = [40 9; 20 5; 60 15];  % [prey0, pred0] rows

figure; hold on;
colors = {'b', 'r', 'g'};
for k = 1:size(initial_conditions, 1)
    prey = zeros(size(t));
    pred = zeros(size(t));
    prey(1) = initial_conditions(k, 1);
    pred(1) = initial_conditions(k, 2);
    for i = 1:length(t)-1
        dprey = a*prey(i) - b*prey(i)*pred(i);
        dpred = -c*pred(i) + d*prey(i)*pred(i);
        prey(i+1) = prey(i) + dprey*dt;
        pred(i+1) = pred(i) + dpred*dt;
    end
    plot(prey, pred, colors{k}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Start: prey=%d, pred=%d', initial_conditions(k,1), initial_conditions(k,2)));
end

xlabel('Prey population');
ylabel('Predator population');
title('Phase Portrait: Predator-Prey Cycles');
legend show;
grid on;






% Exercise 5: Two-species competition model
r1 = 0.5; r2 = 0.5;
K1 = 1000; K2 = 1000;
a12 = 1.2;  % effect of species 2 on species 1
a21 = 0.8;  % effect of species 1 on species 2
dt = 0.01;
t_end = 100;
t = 0:dt:t_end;

N1 = zeros(size(t));
N2 = zeros(size(t));
N1(1) = 50;
N2(1) = 50;

for i = 1:length(t)-1
    dN1 = r1*N1(i)*(1 - (N1(i) + a12*N2(i))/K1);
    dN2 = r2*N2(i)*(1 - (N2(i) + a21*N1(i))/K2);
    N1(i+1) = N1(i) + dN1*dt;
    N2(i+1) = N2(i) + dN2*dt;
end

figure;
plot(t, N1, 'b-', 'LineWidth', 2); hold on;
plot(t, N2, 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Population');
title(sprintf('Competition Model (a_{12}=%.1f, a_{21}=%.1f)', a12, a21));
legend('Species 1', 'Species 2');
grid on;

fprintf('Final: N1 = %.1f, N2 = %.1f\n', N1(end), N2(end));





% Exercise 6a: Logistic map for a few r values
r_values = [2.5, 3.2, 3.5, 3.9];
n_iter = 50;

figure;
for k = 1:length(r_values)
    r = r_values(k);
    N = zeros(1, n_iter);
    N(1) = 0.5;
    for i = 1:n_iter-1
        N(i+1) = r * N(i) * (1 - N(i));
    end
    subplot(2, 2, k);
    plot(N, 'o-', 'LineWidth', 1.5);
    title(sprintf('r = %.1f', r));
    xlabel('Iteration');
    ylabel('N');
    grid on;
end
sgtitle('Logistic Map at Different r Values');

% Exercise 6b: Bifurcation diagram
r_range = 2.5:0.005:4.0;
n_transient = 200;   % discard these iterations (let it settle)
n_keep = 100;        % keep these for plotting

figure; hold on;
for r = r_range
    N = 0.5;
    for i = 1:n_transient
        N = r * N * (1 - N);
    end
    N_vals = zeros(1, n_keep);
    for i = 1:n_keep
        N = r * N * (1 - N);
        N_vals(i) = N;
    end
    plot(r * ones(1, n_keep), N_vals, 'b.', 'MarkerSize', 0.5);
end

xlabel('r');
ylabel('Long-term N values');
title('Bifurcation Diagram of the Logistic Map');
grid on;





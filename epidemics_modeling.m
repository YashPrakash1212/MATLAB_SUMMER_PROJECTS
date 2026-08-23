% Exercise 1: Basic SIR simulation
N = 10000;
I0 = 10;
S0 = N - I0;
R0_init = 0;
beta = 0.3;
gamma = 0.1;
dt = 0.1;
days = 160;
n_steps = round(days/dt);

S = zeros(n_steps+1, 1);
I = zeros(n_steps+1, 1);
R = zeros(n_steps+1, 1);
S(1) = S0; I(1) = I0; R(1) = R0_init;

for t = 1:n_steps
    dS = -beta*S(t)*I(t)/N;
    dI = beta*S(t)*I(t)/N - gamma*I(t);
    dR = gamma*I(t);
    S(t+1) = S(t) + dS*dt;
    I(t+1) = I(t) + dI*dt;
    R(t+1) = R(t) + dR*dt;
end

time = (0:n_steps) * dt;

figure;
plot(time, S, 'b-', 'LineWidth', 2); hold on;
plot(time, I, 'r-', 'LineWidth', 2);
plot(time, R, 'g-', 'LineWidth', 2);
xlabel('Days');
ylabel('People');
title(sprintf('SIR Model (R_0 = %.1f)', beta/gamma));
legend('Susceptible', 'Infected', 'Recovered');
grid on;


% Exercise 2: Vary R0
N = 10000; I0 = 10; gamma = 0.1; dt = 0.1; days = 160;
R0_values = [1.5, 2, 3, 5, 10];
time = (0:round(days/dt)) * dt;

figure; hold on;
for k = 1:length(R0_values)
    beta = R0_values(k) * gamma;
    [S, I, R] = simulate_sir(N, I0, beta, gamma, dt, days);
    plot(time, I, 'LineWidth', 2, 'DisplayName', sprintf('R_0 = %.1f', R0_values(k)));
end
xlabel('Days');
ylabel('Infected');
title('Effect of R_0 on Epidemic Curve');
legend show;
grid on;

function [S, I, R] = simulate_sir(N, I0, beta, gamma, dt, days)
n_steps = round(days/dt);
S = zeros(n_steps+1, 1); I = zeros(n_steps+1, 1); R = zeros(n_steps+1, 1);
S(1) = N - I0; I(1) = I0; R(1) = 0;
for t = 1:n_steps
    dS = -beta*S(t)*I(t)/N;
    dI = beta*S(t)*I(t)/N - gamma*I(t);
    dR = gamma*I(t);
    S(t+1) = S(t) + dS*dt;
    I(t+1) = I(t) + dI*dt;
    R(t+1) = R(t) + dR*dt;
end
end


% Exercise 3: Herd immunity threshold
N = 10000; I0 = 10; beta = 0.3; gamma = 0.1; dt = 0.1; days = 160;
R0_val = beta / gamma;  % = 3
threshold = 1 - 1/R0_val;

vax_rates = linspace(0, 1, 21);
totals = zeros(size(vax_rates));

for k = 1:length(vax_rates)
    totals(k) = simulate_sir_immune(N, I0, beta, gamma, dt, days, vax_rates(k));
end

figure;
plot(vax_rates*100, totals, 'o-', 'Color', [0.5 0 0.5], 'LineWidth', 2); hold on;
xline(threshold*100, 'r--', 'LineWidth', 2, ...
    'Label', sprintf('Herd immunity (%.0f%%)', threshold*100));
xlabel('Vaccination Rate (%)');
ylabel('Total Infected');
title('Vaccination Rate vs Total Infections');
grid on;

function total_new = simulate_sir_immune(N, I0, beta, gamma, dt, days, frac_immune)
n_steps = round(days/dt);
R_init = frac_immune * N;
S_init = N - I0 - R_init;
S = zeros(n_steps+1,1); I = zeros(n_steps+1,1); R = zeros(n_steps+1,1);
S(1) = S_init; I(1) = I0; R(1) = R_init;
for t = 1:n_steps
    dS = -beta*S(t)*I(t)/N;
    dI = beta*S(t)*I(t)/N - gamma*I(t);
    dR = gamma*I(t);
    S(t+1) = S(t) + dS*dt;
    I(t+1) = I(t) + dI*dt;
    R(t+1) = R(t) + dR*dt;
end
total_new = R(end) - R_init;
end



% Exercise 4: Flattening the curve (social distancing)
N = 10000; I0 = 10; gamma = 0.1; dt = 0.1; days = 160;
beta_base = 0.3;
reductions = [0, 0.3, 0.6];
labels = {'No intervention', '30% reduction', '60% reduction'};
time = (0:round(days/dt)) * dt;

figure; hold on;
for k = 1:length(reductions)
    beta = beta_base * (1 - reductions(k));
    [S, I, R] = simulate_sir(N, I0, beta, gamma, dt, days);
    plot(time, I, 'LineWidth', 2, ...
        'DisplayName', sprintf('%s (peak = %.0f)', labels{k}, max(I)));
end
xlabel('Days');
ylabel('Infected');
title('Flattening the Curve: Effect of Social Distancing');
legend show;
grid on;



% Exercise 5: SIRV model with ongoing vaccination
N = 10000; I0 = 10; beta = 0.3; gamma = 0.1; dt = 0.1; days = 160;
time = (0:round(days/dt)) * dt;

[S0v, I0v, R0v] = simulate_sirv(N, I0, beta, gamma, 0, dt, days);      % no vaccination
[S1v, I1v, R1v] = simulate_sirv(N, I0, beta, gamma, 0.02, dt, days);   % 2%/day vaccination

figure;
plot(time, I0v, 'r-', 'LineWidth', 2); hold on;
plot(time, I1v, 'g-', 'LineWidth', 2);
xlabel('Days');
ylabel('Infected');
title('SIRV Model: Impact of a Vaccination Campaign');
legend('No vaccination', 'Vaccination campaign (2%/day)');
grid on;

function [S, I, R] = simulate_sirv(N, I0, beta, gamma, v_rate, dt, days)
n_steps = round(days/dt);
S = zeros(n_steps+1,1); I = zeros(n_steps+1,1); R = zeros(n_steps+1,1);
S(1) = N - I0; I(1) = I0; R(1) = 0;
for t = 1:n_steps
    dS = -beta*S(t)*I(t)/N - v_rate*S(t);
    dI = beta*S(t)*I(t)/N - gamma*I(t);
    dR = gamma*I(t) + v_rate*S(t);
    S(t+1) = S(t) + dS*dt;
    I(t+1) = I(t) + dI*dt;
    R(t+1) = R(t) + dR*dt;
end
end

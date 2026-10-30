
%% DAY 2: ENZYME KINETICS
% Michaelis-Menten and Nonlinear Curve Fitting

clear;
close all;
clc;

rng(42);

%% EXERCISE 1: MICHAELIS-MENTEN CURVE

Vmax_true = 100;
Km_true = 5;

S = linspace(0,50,100);

% Michaelis-Menten equation
mm = @(S,Vmax,Km) (Vmax .* S) ./ (Km + S);

rate = mm(S,Vmax_true,Km_true);

figure;

plot(S,rate,'b','LineWidth',2);
hold on;

yline(Vmax_true,'r--','Vmax');
xline(Km_true,'k--','Km');

xlabel('Substrate Concentration [S]');
ylabel('Reaction Rate');
title('Michaelis-Menten Enzyme Kinetics');

legend('Reaction Rate','Vmax','Km');
grid on;

%% EXERCISE 2: GENERATE NOISY EXPERIMENTAL DATA

S_data = linspace(0.5,50,30);

true_rate = mm(S_data,Vmax_true,Km_true);

% Add experimental noise
noise = 5 * randn(size(S_data));

data_rate = true_rate + noise;

figure;

scatter(S_data,data_rate,40,'filled');
hold on;

plot(S,rate,'r','LineWidth',2);

xlabel('Substrate Concentration [S]');
ylabel('Reaction Rate');

title('Noisy Experimental Enzyme Data');

legend('Experimental Data','True Curve');
grid on;

%% EXERCISE 3: NONLINEAR OPTIMIZATION

% Initial parameter guesses
initial_guess = [80,3];

% Objective function: sum of squared errors
objective = @(p) sum( ...
    (data_rate - mm(S_data,abs(p(1)),abs(p(2)))).^2);

% Optimization using fminsearch
estimated = fminsearch(objective,initial_guess);

Vmax_est = abs(estimated(1));
Km_est = abs(estimated(2));

% Generate fitted curve
S_fit = linspace(0,50,200);

fitted_rate = mm(S_fit,Vmax_est,Km_est);

% Calculate goodness of fit
predicted = mm(S_data,Vmax_est,Km_est);

SSE = sum((data_rate-predicted).^2);

SST = sum((data_rate-mean(data_rate)).^2);

R_squared = 1-SSE/SST;

RMSE = sqrt(mean((data_rate-predicted).^2));

fprintf('\n===== NONLINEAR OPTIMIZATION =====\n');

fprintf('True Vmax: %.4f\n',Vmax_true);
fprintf('Estimated Vmax: %.4f\n',Vmax_est);

fprintf('True Km: %.4f\n',Km_true);
fprintf('Estimated Km: %.4f\n',Km_est);

fprintf('R-squared: %.4f\n',R_squared);
fprintf('RMSE: %.4f\n',RMSE);

% Plot fitted results
figure;

scatter(S_data,data_rate,40,'b','filled');
hold on;

plot(S_fit,fitted_rate,'r','LineWidth',2);

xlabel('Substrate Concentration [S]');
ylabel('Reaction Rate');

title('Nonlinear Optimization: Michaelis-Menten Fit');

legend('Noisy Data','Optimized Curve');
grid on;

saveas(gcf,'michaelis_menten_fit.png');

%% EXERCISE 4: LINEWEAVER-BURK PLOT

% Remove nonpositive noisy reaction rates
valid = data_rate > 0;

x = 1 ./ S_data(valid);
y = 1 ./ data_rate(valid);

% Linear regression
coefficients = polyfit(x,y,1);

slope = coefficients(1);
intercept = coefficients(2);

Vmax_LB = 1/intercept;
Km_LB = slope/intercept;

fprintf('\n===== LINEWEAVER-BURK =====\n');

fprintf('Vmax from Lineweaver-Burk: %.4f\n',Vmax_LB);
fprintf('Km from Lineweaver-Burk: %.4f\n',Km_LB);

% Plot
x_fit = linspace(min(x),max(x),100);

y_fit = polyval(coefficients,x_fit);

figure;

scatter(x,y,40,'b','filled');
hold on;

plot(x_fit,y_fit,'r','LineWidth',2);

xlabel('1 / Substrate Concentration');
ylabel('1 / Reaction Rate');

title('Lineweaver-Burk Plot');

legend('Experimental Data','Linear Fit');
grid on;

saveas(gcf,'lineweaver_burk.png');

%% EXERCISE 5: COMPETITIVE ENZYME INHIBITION

I = 4;
Ki = 2;

% Apparent Km under competitive inhibition
Km_inhibited = Km_true * (1 + I/Ki);

normal_rate = mm(S,Vmax_true,Km_true);

inhibited_rate = mm(S,Vmax_true,Km_inhibited);

figure;

plot(S,normal_rate,'b','LineWidth',2);
hold on;

plot(S,inhibited_rate,'r--','LineWidth',2);

xlabel('Substrate Concentration [S]');
ylabel('Reaction Rate');

title('Competitive Enzyme Inhibition');

legend('Normal Enzyme','Competitive Inhibition');
grid on;

saveas(gcf,'enzyme_inhibition.png');

%% EXERCISE 6: PARAMETER COMPARISON

figure;

bar([Vmax_true,Vmax_est,Km_true,Km_est]);

set(gca,'XTickLabel', ...
    {'True Vmax','Estimated Vmax','True Km','Estimated Km'});

ylabel('Parameter Value');

title('True vs Estimated Enzyme Parameters');

grid on;

fprintf('\n===== PROJECT COMPLETE =====\n');

fprintf('All six enzyme kinetics exercises completed.\n');

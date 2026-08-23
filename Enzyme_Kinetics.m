% Exercise 1: Michaelis-Menten curve
Vmax = 100;
Km = 5;
S = linspace(0, 50, 500);
rate = (Vmax .* S) ./ (Km + S);

figure;
plot(S, rate, 'b-', 'LineWidth', 2); hold on;
yline(Vmax, 'k--', 'LineWidth', 1.5);
xline(Km, 'r--', 'LineWidth', 1.5);
plot(Km, Vmax/2, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(Km + 1, Vmax/2, sprintf('K_m = %d', Km));
text(2, Vmax + 3, sprintf('V_{max} = %d', Vmax));
xlabel('[S] (Substrate concentration)');
ylabel('Reaction rate v');
title('Michaelis-Menten Kinetics');
legend('v vs [S]', 'V_{max}', 'K_m', 'Location', 'southeast');
grid on;
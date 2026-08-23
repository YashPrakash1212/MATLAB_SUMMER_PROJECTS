% Coupled reaction: A -> B -> C
k1 = 0.1;   % rate constant A -> B
k2 = 0.05;  % rate constant B -> C

tspan = [0, 100];       % use square brackets
C0 = [100; 0; 0];       % use square brackets, semicolons stack as a column

[t, C] = ode45(@(t, C) reactionODE(t, C, k1, k2), tspan, C0);

figure;
plot(t, C(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t, C(:,2), 'r-', 'LineWidth', 2);
plot(t, C(:,3), 'g-', 'LineWidth', 2);
xlabel('Time');
ylabel('Concentration');
title('Coupled Reaction A \rightarrow B \rightarrow C');
legend('[A]', '[B]', '[C]');
grid on;

function dCdt = reactionODE(~, C, k1, k2)
A = C(1);
B = C(2);
dA = -k1 * A;
dB = k1 * A - k2 * B;
dC = k2 * B;
dCdt = [dA; dB; dC];   % square brackets here too
end
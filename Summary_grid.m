% Computational Physics: Solving the Unsolvable -- Summary Grid
% One representative plot from each simulation, arranged 2x3:
% Projectile, Pendulum Chaos, and Orbit
% Coupled Springs, Wave, Energy Conservation
%  Each mini-simulation below is a condensed, self-contained version of
%  the physics from the full exercise scripts, just enough to produce
%  one clean plot per panel.

clear; clc; close all;

fig = figure('Name','Computational Physics Summary', 'Position', [80 80 1500 850]);

% Panel 1: Projectile (with air resistance) 
subplot(2,3,1);
g = 9.81; v0 = 50; theta = deg2rad(45); b = 0.1; m = 1; dt = 0.01;
vx = v0*cos(theta); vy = v0*sin(theta); x = 0; y = 0;
xs = x; ys = y;
while y >= 0
    ax = -b*vx/m; ay = -g - b*vy/m;
    vx = vx + ax*dt; vy = vy + ay*dt;
    x = x + vx*dt; y = y + vy*dt;
    xs(end+1) = x; ys(end+1) = y; 
end
plot(xs, ys, 'b-', 'LineWidth', 2); grid on;
xlabel('x (m)'); ylabel('y (m)'); title('Projectile with Air Resistance');

% Panel 2: Pendulum Chaos (double pendulum)
subplot(2,3,2);
m1 = 1; m2 = 1; L1 = 1; L2 = 1; gP = 9.8; dtP = 0.0005; TmaxP = 15;
NP = round(TmaxP/dtP);
th1 = deg2rad(120); th2 = deg2rad(100); w1 = 0; w2 = 0;
x2h = zeros(NP+1,1); y2h = zeros(NP+1,1);
x1c = L1*sin(th1); y1c = -L1*cos(th1);
x2h(1) = x1c + L2*sin(th2); y2h(1) = y1c - L2*cos(th2);
for i = 1:NP
    denom = 2*m1 + m2 - m2*cos(2*th1 - 2*th2);
    num1 = -gP*(2*m1+m2)*sin(th1) - m2*gP*sin(th1-2*th2) ...
           - 2*sin(th1-th2)*m2*(w2^2*L2 + w1^2*L1*cos(th1-th2));
    a1 = num1/(L1*denom);
    num2 = 2*sin(th1-th2)*(w1^2*L1*(m1+m2) + gP*(m1+m2)*cos(th1) + w2^2*L2*m2*cos(th1-th2));
    a2 = num2/(L2*denom);
    w1 = w1 + a1*dtP; w2 = w2 + a2*dtP;
    th1 = th1 + w1*dtP; th2 = th2 + w2*dtP;
    x1c = L1*sin(th1); y1c = -L1*cos(th1);
    x2h(i+1) = x1c + L2*sin(th2); y2h(i+1) = y1c - L2*cos(th2);
end
plot(x2h, y2h, 'm-', 'LineWidth', 0.5); axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)'); title('Double Pendulum Chaos');

% Panel 3: Orbit (elliptical) 
subplot(2,3,3);
G = 6.674e-11; M = 1e30; GM = G*M;
r0 = [1e11, 0]; v0o = [0, 2e4]; dtO = 1000;
rV = norm(r0); vV = norm(v0o);
a_semi = 1/(2/rV - vV^2/GM);
T_orbit = 2*pi*sqrt(a_semi^3/GM);
NO = round(1.2*T_orbit/dtO);
rx = zeros(NO+1,1); ry = zeros(NO+1,1); vx_o = zeros(NO+1,1); vy_o = zeros(NO+1,1);
rx(1) = r0(1); ry(1) = r0(2); vx_o(1) = v0o(1); vy_o(1) = v0o(2);
for i = 1:NO
    rMag = sqrt(rx(i)^2 + ry(i)^2);
    axo = -GM*rx(i)/rMag^3; ayo = -GM*ry(i)/rMag^3;
    vx_o(i+1) = vx_o(i) + axo*dtO; vy_o(i+1) = vy_o(i) + ayo*dtO;
    rx(i+1) = rx(i) + vx_o(i+1)*dtO; ry(i+1) = ry(i) + vy_o(i+1)*dtO;
end
plot(rx, ry, 'b-', 'LineWidth', 1.5); hold on;
plot(0, 0, 'yo', 'MarkerSize', 10, 'MarkerFaceColor', 'y'); hold off;
axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)'); title('Elliptical Orbit');

% Panel 4: Coupled Springs (normal-mode beating) 
subplot(2,3,4);
mS = 1; k = 1; k12 = 0.15; dtS = 0.01; TmaxS = 100;
NS = round(TmaxS/dtS);
state = [1; 0; 0; 0];   % [x1, x2, v1, v2] -- pluck mass 1 only
X1 = zeros(NS+1,1); X2 = zeros(NS+1,1); tS = zeros(NS+1,1);
X1(1) = state(1); X2(1) = state(2);
springDeriv = @(s) [s(3); s(4); ...
    (-k*s(1) - k12*(s(1)-s(2)))/mS; ...
    (-k*s(2) - k12*(s(2)-s(1)))/mS];
for i = 1:NS
    k1v = springDeriv(state);
    k2v = springDeriv(state + dtS/2*k1v);
    k3v = springDeriv(state + dtS/2*k2v);
    k4v = springDeriv(state + dtS*k3v);
    state = state + (dtS/6)*(k1v + 2*k2v + 2*k3v + k4v);
    X1(i+1) = state(1); X2(i+1) = state(2); tS(i+1) = tS(i) + dtS;
end
plot(tS, X1, 'b-', 'LineWidth', 1.2); hold on;
plot(tS, X2, 'r-', 'LineWidth', 1.2); hold off; grid on;
xlabel('Time (s)'); ylabel('Displacement (m)'); title('Coupled Springs: Energy Beating');
legend('Mass 1', 'Mass 2', 'Location', 'best', 'FontSize', 7);

% Panel 5: Wave on a String 
subplot(2,3,5);
Nx = 100; L = 1; xw = linspace(0, L, Nx); dxw = xw(2)-xw(1);
cW = 1; CourantW = 0.5; dtw = CourantW*dxw/cW;
A = 0.1;
y0w = zeros(1,Nx);
for i = 1:Nx
    if xw(i) <= L/2, y0w(i) = 2*A*xw(i)/L; else, y0w(i) = 2*A*(L-xw(i))/L; end
end
Nsteps = round(0.6/dtw);
Y = zeros(Nx, Nsteps+1); Y(:,1) = y0w(:);
Cc = cW*dtw/dxw;
Y(2:Nx-1,2) = Y(2:Nx-1,1) + 0.5*Cc^2*(Y(3:Nx,1)-2*Y(2:Nx-1,1)+Y(1:Nx-2,1));
for n = 2:Nsteps
    Y(2:Nx-1,n+1) = 2*Y(2:Nx-1,n) - Y(2:Nx-1,n-1) + Cc^2*(Y(3:Nx,n)-2*Y(2:Nx-1,n)+Y(1:Nx-2,n));
end
plot(xw, Y(:,1), 'b-', 'LineWidth', 1.5); hold on;
plot(xw, Y(:,round(Nsteps*0.4)), 'g-', 'LineWidth', 1.5);
plot(xw, Y(:,end), 'r-', 'LineWidth', 1.5); hold off;
grid on; ylim([-0.13 0.13]);
xlabel('x (m)'); ylabel('y (m)'); title('Wave on a String');
legend('t=0','later','later still', 'Location', 'best', 'FontSize', 7);

% Panel 6: Energy Conservation (Euler vs RK4) 
subplot(2,3,6);
r0e = [1e11, 0]; v0e = [0, 3e4]; dte = 1000;
rVe = norm(r0e); vVe = norm(v0e);
a_e = 1/(2/rVe - vVe^2/GM);
Te = 2*pi*sqrt(a_e^3/GM);
Ne = round(3*Te/dte);

% Naive Euler
rxE=zeros(Ne+1,1); ryE=zeros(Ne+1,1); vxE=zeros(Ne+1,1); vyE=zeros(Ne+1,1);
rxE(1)=r0e(1); ryE(1)=r0e(2); vxE(1)=v0e(1); vyE(1)=v0e(2);
for i=1:Ne
    rMag=sqrt(rxE(i)^2+ryE(i)^2);
    axE=-GM*rxE(i)/rMag^3; ayE=-GM*ryE(i)/rMag^3;
    rxE(i+1)=rxE(i)+vxE(i)*dte; ryE(i+1)=ryE(i)+vyE(i)*dte;
    vxE(i+1)=vxE(i)+axE*dte;   vyE(i+1)=vyE(i)+ayE*dte;
end

% RK4
stateR=[r0e(1); r0e(2); v0e(1); v0e(2)];
rxR=zeros(Ne+1,1); ryR=zeros(Ne+1,1); vxR=zeros(Ne+1,1); vyR=zeros(Ne+1,1);
rxR(1)=stateR(1); ryR(1)=stateR(2); vxR(1)=stateR(3); vyR(1)=stateR(4);
deriv = @(s) [s(3); s(4); -GM*s(1)/norm(s(1:2))^3; -GM*s(2)/norm(s(1:2))^3];
for i=1:Ne
    k1o=deriv(stateR); k2o=deriv(stateR+dte/2*k1o);
    k3o=deriv(stateR+dte/2*k2o); k4o=deriv(stateR+dte*k3o);
    stateR = stateR + (dte/6)*(k1o+2*k2o+2*k3o+k4o);
    rxR(i+1)=stateR(1); ryR(i+1)=stateR(2); vxR(i+1)=stateR(3); vyR(i+1)=stateR(4);
end

tvec = (0:Ne)*dte/86400;
E_E = 0.5*(vxE.^2+vyE.^2) - GM./sqrt(rxE.^2+ryE.^2);
E_R = 0.5*(vxR.^2+vyR.^2) - GM./sqrt(rxR.^2+ryR.^2);

plot(tvec, E_E, 'r-', 'LineWidth', 1.3); hold on;
plot(tvec, E_R, 'b-', 'LineWidth', 1.3); hold off; grid on;
xlabel('Time (days)'); ylabel('Specific Energy (J/kg)');
title('Energy Conservation: Euler vs RK4');
legend('Euler (drifts)', 'RK4 (stays flat)', 'Location', 'best', 'FontSize', 7);

% Overall title 
sgtitle('Computational Physics: Solving the Unsolvable', 'FontSize', 18, 'FontWeight', 'bold');


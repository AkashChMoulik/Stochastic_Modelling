% MATLAB Program to simulate SDOF Duffing-Van der Pol Oscillator

% Define parameters
mu = 1.0;         % Coefficient of the nonlinear damping term
omega0 = 1.0;     % Natural frequency of the linear part of the system
alpha = 1.0;      % Coefficient of the nonlinear stiffness term
F = 0.5;          % Amplitude of the external forcing
omega = 1.2;      % Frequency of the external forcing

% Time span for simulation
tspan = [0 100];

% Initial conditions [x0, v0]
initial_conditions = [0.1, 0]; % Initial displacement and velocity

% Define the Duffing-Van der Pol ODE as a function
duffing_vdp = @(t, y) [y(2); 
                       F*cos(omega*t) - mu*(y(1)^2 - 1)*y(2) - omega0^2*y(1) - alpha*y(1)^3];

% Solve the ODE using ode45 solver
[t, Y] = ode45(duffing_vdp, tspan, initial_conditions);

% Extract displacement and velocity
x = Y(:, 1);
v = Y(:, 2);

% Plot the results
figure;
subplot(2, 1, 1);
plot(t, x, 'b');
xlabel('Time');
ylabel('Displacement');
title('Displacement vs Time');
grid on;

subplot(2, 1, 2);
plot(x, v, 'r');
xlabel('Displacement');
ylabel('Velocity');
title('Phase Portrait');
grid on;

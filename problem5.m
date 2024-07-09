% MATLAB Program to simulate SDOF Duffing Oscillator

% Define parameters
zeta = 0.05;        % Damping ratio
omega0 = 1.0;       % Natural frequency of the linear part of the system
alpha = 1.0;        % Coefficient of the nonlinear stiffness term
F = 0.5;            % Amplitude of the external forcing
omega = 1.2;        % Frequency of the external forcing

% Time span for simulation
tspan = [0 100];

% Initial conditions [x0, v0]
initial_conditions = [0.1, 0]; % Initial displacement and velocity

% Define the Duffing ODE as a function
duffing = @(t, y) [y(2); 
                   F*cos(omega*t) - 2*zeta*omega0*y(2) - omega0^2*y(1) - alpha*y(1)^3];

% Solve the ODE using ode45 solver
[t, Y] = ode45(duffing, tspan, initial_conditions);

% Extract displacement and velocity
x = Y(:, 1);
v = Y(:, 2);

% Plot the results
figure;

% Displacement vs Time
subplot(2, 1, 1);
plot(t, x, 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Displacement');
title('Displacement vs Time');
grid on;

% Phase Portrait
subplot(2, 1, 2);
plot(x, v, 'r', 'LineWidth', 1.5);
xlabel('Displacement');
ylabel('Velocity');
title('Phase Portrait');
grid on;

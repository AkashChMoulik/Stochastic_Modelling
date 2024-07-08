% Parameters
lambda = 2;  % Drift coefficient
mu = 1;      % Diffusion coefficient
X0 = 1;      % Initial value
T = 1;       % Final time
N = 1000;    % Number of time steps
dt = T / N;  % Time step size
t = linspace(0, T, N + 1);  % Time vector

% Brownian motion
W = zeros(1, N + 1);  % Initialize Brownian motion
W(2:end) = cumsum(sqrt(dt) * randn(1, N));  % Generate Brownian increments

% Exact solution
X_exact = X0 * exp((lambda - 0.5 * mu^2) * t + mu * W);  % Exact solution

% Euler's method
X_euler = zeros(1, N + 1);  % Initialize Euler approximation
X_euler(1) = X0;
for i = 1:N
    X_euler(i + 1) = X_euler(i) + lambda * X_euler(i) * dt + mu * X_euler(i) * (W(i + 1) - W(i));
end

% Plot results
figure;
plot(t, X_exact, 'r', 'DisplayName', 'Exact Solution');
hold on;
plot(t, X_euler, 'b--', 'DisplayName', 'Euler Approximation');
xlabel('Time');
ylabel('X(t)');
legend;
title('Comparison of Exact Solution and Euler Approximation');

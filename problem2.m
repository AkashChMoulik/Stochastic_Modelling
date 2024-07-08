% Parameters
X0 = 1;
T = 1;
N = 1000;
dt = T/N;
num_simulations = 100;
lambda_range = [1.8, 2.1];
mu = 1;

% Time vector
t = linspace(0, T, N+1);

% Exact solution
lambda_exact = 2;
W = [0 cumsum(sqrt(dt)*randn(1, N))]; % Brownian motion
X_exact = X0 * exp((lambda_exact - mu^2 / 2) * t + mu * W);

% Monte Carlo simulation
X_ensemble = zeros(num_simulations, N+1);
for j = 1:num_simulations
    lambda = rand() * (lambda_range(2) - lambda_range(1)) + lambda_range(1);
    W = [0 cumsum(sqrt(dt)*randn(1, N))]; % Brownian motion for each simulation
    X = zeros(1, N+1);
    X(1) = X0;
    for i = 1:N
        X(i+1) = X(i) + lambda * X(i) * dt + mu * X(i) * (W(i+1) - W(i));
    end
    X_ensemble(j, :) = X;
end

% Ensemble average
X_avg = mean(X_ensemble);

% Plot results
figure;
plot(t, X_exact, 'r', 'DisplayName', 'Exact Solution');
hold on;
plot(t, X_avg, 'b--', 'DisplayName', 'Ensemble Average');
xlabel('Time');
ylabel('X(t)');
legend;
title('Comparison of Exact Solution and Ensemble Average (Uniform \lambda)');

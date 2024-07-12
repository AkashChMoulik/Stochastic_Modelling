
% Parameters
mu = 1.0;
sigma = 0.1;
T = 10; % Total time
num_steps = 10^6; % Number of time steps

% Time step sizes for each method
dt_EulerMaruyama = 1e-6;
dt_Milstein = 1e-4;
dt_Taylor15 = 1e-2;

% Initial conditions
X0 = 1.0;
Y0 = 1.0;

% Preallocate arrays for storing results
X_EM = zeros(num_steps+1, 1);
Y_EM = zeros(num_steps+1, 1);
X_Milstein = zeros(num_steps+1, 1);
Y_Milstein = zeros(num_steps+1, 1);
X_Taylor15 = zeros(num_steps+1, 1);
Y_Taylor15 = zeros(num_steps+1, 1);

% Initial conditions
X_EM(1) = X0;
Y_EM(1) = Y0;
X_Milstein(1) = X0;
Y_Milstein(1) = Y0;
X_Taylor15(1) = X0;
Y_Taylor15(1) = Y0;

% Wiener increments
dW_EM = sqrt(dt_EulerMaruyama) * randn(num_steps, 1);
dW_Milstein = sqrt(dt_Milstein) * randn(num_steps, 1);
dW_Taylor15 = sqrt(dt_Taylor15) * randn(num_steps, 1);

% Simulation loop
for i = 1:num_steps
    % Euler-Maruyama method
    X_EM(i+1) = X_EM(i) + Y_EM(i) * dt_EulerMaruyama;
    Y_EM(i+1) = Y_EM(i) + (mu * (1 - X_EM(i)^2) * Y_EM(i) - X_EM(i)) * dt_EulerMaruyama + sigma * dW_EM(i);
    
    % Milstein method
    X_Milstein(i+1) = X_Milstein(i) + Y_Milstein(i) * dt_Milstein;
    Y_Milstein(i+1) = Y_Milstein(i) + (mu * (1 - X_Milstein(i)^2) * Y_Milstein(i) - X_Milstein(i)) * dt_Milstein + ...
                      0.5 * sigma * (Y_Milstein(i) * (2 * X_Milstein(i) * dW_Milstein(i) + sigma * dW_Milstein(i)^2) ...
                      + (mu * (1 - X_Milstein(i)^2) - 1) * dt_Milstein);
    
    % Taylor 1.5 method
    X_Taylor15(i+1) = X_Taylor15(i) + Y_Taylor15(i) * dt_Taylor15;
    Y_Taylor15(i+1) = Y_Taylor15(i) + (mu * (1 - X_Taylor15(i)^2) * Y_Taylor15(i) - X_Taylor15(i)) * dt_Taylor15 + ...
                      0.5 * sigma * (Y_Taylor15(i) * (2 * X_Taylor15(i) * dW_Taylor15(i) + sigma * dW_Taylor15(i)^2) ...
                      + (mu * (1 - X_Taylor15(i)^2) - 1) * dt_Taylor15);
end

% Plot results
figure;
subplot(2, 1, 1);
plot(linspace(0, T, num_steps+1), X_EM, 'b', 'LineWidth', 1.5);
hold on;
plot(linspace(0, T, num_steps+1), X_Milstein, 'r', 'LineWidth', 1.5);
plot(linspace(0, T, num_steps+1), X_Taylor15, 'g', 'LineWidth', 1.5);
xlabel('Time');
ylabel('X(t)');
title('Comparison of X(t) using Different Methods');
legend('Euler-Maruyama', 'Milstein', 'Taylor 1.5');

subplot(2, 1, 2);
plot(linspace(0, T, num_steps+1), Y_EM, 'b', 'LineWidth', 1.5);
hold on;
plot(linspace(0, T, num_steps+1), Y_Milstein, 'r', 'LineWidth', 1.5);
plot(linspace(0, T, num_steps+1), Y_Taylor15, 'g', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Y(t)');
title('Comparison of Y(t) using Different Methods');
legend('Euler-Maruyama', 'Milstein', 'Taylor 1.5');
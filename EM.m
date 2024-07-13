% Parameters
mu = 1.0;
sigma = 0.1;
T = 1; % Total time changed to 1
num_steps = 10^6; % Number of time steps

% Time step sizes for each method
dt_EulerMaruyama = T / num_steps;
dt_Milstein = T / num_steps;
dt_Taylor15 = T / num_steps;

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
                      sigma * dW_Milstein(i) + 0.5 * sigma^2 * dW_Milstein(i) * (dW_Milstein(i) - dt_Milstein);
    
    % Taylor 1.5 method
    X_Taylor15(i+1) = X_Taylor15(i) + Y_Taylor15(i) * dt_Taylor15;
    Y_Taylor15(i+1) = Y_Taylor15(i) + (mu * (1 - X_Taylor15(i)^2) * Y_Taylor15(i) - X_Taylor15(i)) * dt_Taylor15 + ...
                      sigma * dW_Taylor15(i) + 0.5 * sigma^2 * dW_Taylor15(i) * (dW_Taylor15(i) - dt_Taylor15);
end

% Plot results
time = linspace(0, T, num_steps+1);

figure;
subplot(2, 1, 1);
plot(time, X_EM, 'b', 'LineWidth', 2);
hold on;
plot(time, X_Milstein, 'r', 'LineWidth', 2);
plot(time, X_Taylor15, 'g', 'LineWidth', 2);
xlim([0 1]);
xlabel('Time');
ylabel('X(t)');
title('Comparison of X(t) using Different Methods');
legend('Euler-Maruyama', 'Milstein', 'Taylor 1.5');

subplot(2, 1, 2);
plot(time, Y_EM, 'b', 'LineWidth', 2);
hold on;
plot(time, Y_Milstein, 'r', 'LineWidth', 2);
plot(time, Y_Taylor15, 'g', 'LineWidth', 2);
xlim([0 1]);
xlabel('Time');
ylabel('Y(t)');
title('Comparison of Y(t) using Different Methods');
legend('Euler-Maruyama', 'Milstein', 'Taylor 1.5');

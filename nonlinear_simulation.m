% MATLAB code for numerical simulation of nonlinear structural systems

% Define parameters for Duffing-Van der Pol oscillator
alpha = 1;
beta = -1;
delta = 0.5;
gamma = 0.5;
omega = 1.2;

% Time parameters
T = 100; % Total time
dt = 0.01; % Time step
N = T/dt; % Number of time steps

% Initial conditions
x0 = 0.1;
v0 = 0.1;

% Array of sigma values
sigma_values = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35];

% Loop through different sigma values
for sigma_idx = 1:length(sigma_values)
    sigma = sigma_values(sigma_idx);

    % Initialize arrays to store results
    t = (0:dt:T-dt)';
    x = zeros(N, 1);
    v = zeros(N, 1);

    % Set initial conditions
    x(1) = x0;
    v(1) = v0;

    % Integrate using Weak Taylor 3.0 method
    for i = 1:N-1
        W = sqrt(dt) * randn; % Wiener process increment
        x(i+1) = x(i) + v(i)*dt;
        v(i+1) = v(i) + (-delta*v(i) - alpha*x(i) - beta*x(i)^3 + gamma*cos(omega*t(i)))*dt + sigma*W;
    end

    % Plot results
    figure;
    subplot(2,1,1);
    plot(t, x);
    title(['Displacement of SDOF Duffing-Van der Pol Oscillator (\sigma = ', num2str(sigma), ')']);
    xlabel('Time');
    ylabel('Displacement');

    subplot(2,1,2);
    plot(t, v);
    title(['Velocity of SDOF Duffing-Van der Pol Oscillator (\sigma = ', num2str(sigma), ')']);
    xlabel('Time');
    ylabel('Velocity');
end

% Parameters for Two-DOF system with coupled Duffing-Van der Pol oscillator
m1 = 1; m2 = 1;
c1 = 0.1; c2 = 0.1;
k1 = 1; k2 = 1;
alpha2 = 1; beta2 = -1;

% Initial conditions
x1_0 = 0.1; v1_0 = 0.1;
x2_0 = 0.1; v2_0 = 0.1;

% Loop through different sigma values for the Two-DOF system
for sigma_idx = 1:length(sigma_values)
    sigma = sigma_values(sigma_idx);

    % Initialize arrays to store results
    t = (0:dt:T-dt)';
    x1 = zeros(N, 1);
    v1 = zeros(N, 1);
    x2 = zeros(N, 1);
    v2 = zeros(N, 1);

    % Set initial conditions
    x1(1) = x1_0;
    v1(1) = v1_0;
    x2(1) = x2_0;
    v2(1) = v2_0;

    % Integrate using Weak Taylor 3.0 method
    for i = 1:N-1
        W1 = sqrt(dt) * randn; % Wiener process increment for first DOF
        W2 = sqrt(dt) * randn; % Wiener process increment for second DOF

        x1(i+1) = x1(i) + v1(i)*dt;
        v1(i+1) = v1(i) + (-c1*v1(i) - k1*x1(i) - alpha*x1(i)^3 + gamma*cos(omega*t(i)))*dt + sigma*W1;

        x2(i+1) = x2(i) + v2(i)*dt;
        v2(i+1) = v2(i) + (-c2*v2(i) - k2*x2(i) - alpha2*x2(i)^3 + gamma*cos(omega*t(i)))*dt + sigma*W2;
    end

    % Plot results for Two-DOF system
    figure;
    subplot(2,1,1);
    plot(t, x1, t, x2);
    title(['Displacement of Two-DOF Nonlinear System (\sigma = ', num2str(sigma), ')']);
    xlabel('Time');
    ylabel('Displacement');
    legend('x1', 'x2');

    subplot(2,1,2);
    plot(t, v1, t, v2);
    title(['Velocity of Two-DOF Nonlinear System (\sigma = ', num2str(sigma), ')']);
    xlabel('Time');
    ylabel('Velocity');
    legend('v1', 'v2');
end

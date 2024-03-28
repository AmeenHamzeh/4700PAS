function [] = MC()
    % Constants
    q_e = 1.60217653e-19;     % electron charge
    m_e = 9.10938215e-31;     % electron mass
    dt = 1e-14;               % time step size
    num_steps = 1000;         % number of time steps
    N = 10;                  % number of electrons
    scattering_probability = 0.05;
    force = 1e-12;            % constant force for simplicity

    % Initialize positions and velocities
    x = zeros(N, num_steps);
    v = zeros(N, num_steps);

    % Time evolution
    for n = 1:N
        for t = 2:num_steps
            if rand() < scattering_probability
                v(n, t-1) = -0.25*v(n, t-1); % scattering event resets the velocity
            end
            v(n, t) = v(n, t-1) + (force / m_e)*dt;
            x(n, t) = x(n, t-1) + v(n, t) * dt;
        end
    end

    % Plotting
    figure;
    subplot(2,1,1); % Position plot
    plot(0:dt:(num_steps-1)*dt, x);
    xlabel('Time (s)');
    ylabel('Position (m)');
    title('Position vs. Time');

    subplot(2,1,2); % Velocity plot
    plot(0:dt:(num_steps-1)*dt, v);
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    title('Velocity vs. Time');

    % Drift velocity
    drift_velocity = mean(v(:, end));
    sgtitle(['Drift Velocity: ', num2str(drift_velocity), ' m/s']);

    % Make the simulation into a movie if desired
    % (code for creating a movie would go here)
end

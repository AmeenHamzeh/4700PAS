% Grid parameters
nx = 100;    % Number of grid points along x-axis
ny = 100;    % Number of grid points along y-axis
Lx = 1;     % Length of the domain along x-axis
Ly = 1;     % Length of the domain along y-axis

dx = Lx / (nx - 1); % Grid spacing along x-axis
dy = Ly / (ny - 1); % Grid spacing along y-axis

% Boundary conditions
V = zeros(nx, ny);   % Initialize potential matrix with boundary conditions
V(end, :) = 1;      % Right boundary condition
V(1,:) = 1;        % Left boundary condition
V(:, end) = 0;      % Top boundary condition
V(:,1) = 0;  % Bottom boundary condition

% Jacobi Iteration
max_iter = 10000;   % Maximum number of iterations
tolerance = 1e-6;   % Tolerance for convergence
err = Inf;          % Initialize error

for iter = 1:max_iter
    V_old = V;
    % for i = 2:ny-1
    %     % for j = 2:nx-1
    %     %     % Laplace's equation: V_xx + V_yy = 0
    %     %     V(i, j) = 0.25 * (V_old(i+1, j) + V_old(i-1, j) + V_old(i, j+1) + V_old(i, j-1));
    %     % end
    % end
    V = imboxfilt(V,3);
   % Initialize potential matrix with boundary conditions
    V(end, :) = 1;      % Right boundary condition
    V(1,:) = 1;        % Left boundary condition
    V(:, end) = 0;      % Top boundary condition
    V(:,1) = 0;  % Bottom boundary condition
    %V(end,:) = V(end-1,:);
    %V(1,:) = V(2,:);
    if mod(iter,100)==0

        
        surf(V');
        view(45, 135); % Rotate by 90 degrees on z-axis
        pause(0.05)
    end
    % Check for convergence
    err = max(abs(V(:) - V_old(:)));
    if err < tolerance
        break;
    end
end

% Plot the potential field
[X, Y] = meshgrid(0:dx:Lx, 0:dy:Ly);
[Ex, Ey] = gradient(V);
quiver(-Ey',-Ex', 20);
figure
contourf(X, Y, V', 20);
view(180, 90);
figure
surf(V');
view(45, 135);
colorbar;
xlabel('x');
ylabel('y');
title('2D Laplace Solver using Jacobi Iteration');
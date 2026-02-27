% RP Activity 1
% Author: Joel Minj
% Steady-state two-group neutron diffusion equations solver

clear; clc; close all;

% 1. Physical Parameters

Lx = 40;              % width in x (cm)
Ly = 40;              % height in y (cm)
S_strength = 5e5;     % Fast neutron source at the center (neutrons/cm^3/s)

% Group 1 (Fast Neutrons) Parameters
D1 = 1.5;             % Fast diffusion coefficient (cm)
Sigma_a1 = 0.010;     % Fast absorption cross-section (cm^-1)
Sigma_s12 = 0.025;    % Down-scattering cross-section (Fast -> Thermal) (cm^-1)
Sigma_R1 = Sigma_a1 + Sigma_s12; % Total removal from Group 1

% Group 2 (Thermal Neutrons) Parameters
D2 = 0.4;             % Thermal diffusion coefficient (cm)
Sigma_a2 = 0.100;     % Thermal absorption cross-section (cm^-1)


% 2. Spatial Grid Setup (Solving for inside nodes only)

% Assumption: Dirichlet boundary conditions (flux = 0 at the edges)
% Therefore, we only solve for the unknown interior points.
Nx = 200;             % Number of interior nodes in X
Ny = 200;             % Number of interior nodes in Y
N = Nx * Ny;          % Total number of unknown variables per group

dx = Lx / (Nx + 1);   % Mesh spacing in X
dy = Ly / (Ny + 1);   % Mesh spacing in Y

% Spatial coordinates for plotting (including boundaries)
x = linspace(0, Lx, Nx + 2);
y = linspace(0, Ly, Ny + 2);
[X, Y] = meshgrid(x, y);

% 3. Matrix Assembly using Kronecker Tensor Product (kron)

% Create 1D 2nd-Derivative Matrices (Tridiagonal)
e_x = ones(Nx, 1);
e_y = ones(Ny, 1);
D_xx = spdiags([-e_x, 2*e_x, -e_x], [-1, 0, 1], Nx, Nx) / dx^2;
D_yy = spdiags([-e_y, 2*e_y, -e_y], [-1, 0, 1], Ny, Ny) / dy^2;

% Identity matrices for the orthogonal directions
I_x = speye(Nx);
I_y = speye(Ny);

% 2D Laplacian Matrix (Pentadiagonal)
% grad^2 = d^2/dx^2 + d^2/dy^2
Laplacian_2D = kron(I_y, D_xx) + kron(D_yy, I_x); 

% System Matrices for Group 1 and Group 2
% Equation: -D * Laplacian + Sigma * Identity
I_2D = speye(N);
A1 = D1 * Laplacian_2D + Sigma_R1 * I_2D;
A2 = D2 * Laplacian_2D + Sigma_a2 * I_2D;

% 4. Source Definition
% The source only emits fast neutrons (center).

S1 = zeros(N, 1);

% Find the 1D index of the 2D center coordinate
center_i = round(Nx / 2);
center_j = round(Ny / 2);
center_idx = center_i + (center_j - 1) * Nx; 

% Distribute the source strength volumetrically
S1(center_idx) = S_strength / (dx * dy);

% 5. Solve the Coupled System (Sequential Solve)

% Step A: Solve Fast Group first (since its source is independent)
phi1_vec = A1 \ S1;

% Step B: Calculate Thermal Source
% Thermal neutrons are born entirely from fast neutrons slowing down
S2 = Sigma_s12 * phi1_vec;

% Step C: Solve Thermal Group
phi2_vec = A2 \ S2;


% 6. Plotting

% Reshape the 1D solution vectors back into 2D interior grids (Nx by Ny)
phi1_interior = reshape(phi1_vec, Nx, Ny)';
phi2_interior = reshape(phi2_vec, Nx, Ny)';

% Pad with zeros to account for the boundaries (Dirichlet BCs)
phi1_2D = zeros(Nx + 2, Ny + 2);
phi1_2D(2:end-1, 2:end-1) = phi1_interior;

phi2_2D = zeros(Nx + 2, Ny + 2);
phi2_2D(2:end-1, 2:end-1) = phi2_interior;

% Plot
figure('Position', [100, 100, 1400, 600]);
set(gcf, 'DefaultAxesFontSize', 11);

% Fast Flux Heatmap
subplot(1, 2, 1);
contourf(X, Y, phi1_2D, 50, 'LineColor', 'none');
colormap jet; 
cbar1 = colorbar;
ylabel(cbar1, 'Flux [$n/cm^2/s$]', 'FontSize', 11, 'Interpreter', 'latex');
title('Fast Neutron Flux ($\phi_1$)', 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold');
xlabel('$x$ [$cm$]', 'FontSize', 12, 'Interpreter', 'latex'); 
ylabel('$y$ [$cm$]', 'FontSize', 12, 'Interpreter', 'latex');
axis equal tight;
set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');

% Thermal Flux Heatmap
subplot(1, 2, 2);
contourf(X, Y, phi2_2D, 50, 'LineColor', 'none');
colormap jet;
cbar2 = colorbar;
ylabel(cbar2, 'Flux [$n/cm^2/s$]', 'FontSize', 11, 'Interpreter', 'latex');
title('Thermal Neutron Flux ($\phi_2$)', 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold');
xlabel('$x$ [$cm$]', 'FontSize', 12, 'Interpreter', 'latex'); 
ylabel('$y$ [$cm$]', 'FontSize', 12, 'Interpreter', 'latex');
axis equal tight;
set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');

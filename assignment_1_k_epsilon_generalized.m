% Copy Left! by Vagesh D. Narasimhamurthy
% DNS data at Re_delta=7890, Re_tau=395 (Moser, Kim & Mansour, PoF, 1999).
% All quantites are normalized by u_tau and nu unless stated otherwise.
% Delta denotes the channel half-width.

% Author: Nirupam Pal
% Date: 18/03/2024
% Description: This program is used for solving the fully developed channel
%              flow using k-ε model

clc; close all; clear;
image_type = 'epsc';

% Open files for data loggging
file = fopen("solution_k_epsilon.txt", "w");
file1 = fopen('simulated_data_k_epsilon.dat', 'w');

fprintf('*-----------------------------------------*-----------------------------------------*\n');
fprintf('|                                 Author: Nirupam Pal                               |\n');
fprintf('|                                   Date: 18/03/2024                                |\n');
fprintf('|     Description: This program is used for solving the fully developed channel     |\n');
fprintf('|                                 flow using k-ε model                              |\n');
fprintf('*-----------------------------------------*---------------------------------------- *\n\n');

% Read DNS data [half-channel is given (till centerline)]
load y_dns.dat
load u_dns.dat
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat
load uv_dns.dat
load dns_data.dat

% Flow Parameters
nu = 1/395;
ustar = 1;
rho = 1;

% k-ε model constants
c_mu = 0.09;
c1_eps = 1.44;
c2_eps = 1.92;
sigma_k = 1;
sigma_eps = 1.3;
kappa = 0.41; % Von-Karman constant

% Residual error limit
residue_limit = 10^(-4);

% Under relaxation factor
urf = 0.8;

%Grid (based on DNS data)
ny = length(y_dns); % No. of grid points
dy = zeros(ny-1, 1); % Array for increments
% Calutating 'dy' using the DNS data
for i = 1:ny-1
    dy(i) = y_dns(i+1) - y_dns(i);
end

% Creating the control surfaces (Method-1)
y_surface1 = zeros(ny-1, 1);
y_surface1(ny-1) = y_dns(ny);
y_surface1(1) = y_dns(1);
for i=2:ny-2
    y_surface1(i) = y_dns(i) + dy(i)/2;
end

% Creating the control surfaces (Method-2)
% Method-2 is used in the program (Both works well)
y_surface = zeros(ny-1, 1);
y_surface(ny-1) = y_dns(ny);
y_surface(1) = y_dns(1);
for i=2:ny-2
    y_surface(i) = y_surface(i-1) + (y_dns(i)-y_surface(i-1))*2;
end

%Initial conditions for old & new variables (U, dUdy, k, eps, nu_t, residue, ...)
U = zeros(ny, 1);
dUdy = ones(ny-1, 1);
k = 10*ones(ny, 1);
eps = 10*ones(ny, 1);
nu_t = 2*ones(ny, 1);
Pk = zeros(ny, 1);
y_plus = zeros(ny, 1);
U_res = 2*residue_limit; % Initial U residual
k_res = 2*residue_limit; % Initial k residual
eps_res = 2*residue_limit; % Initial ε residual
res = [U_res; k_res; eps_res];

% Boundary conditions
U(1) = 0;
nu_t(1) = 0;
k(1) = 0;
Pk(1) = 0;
dUdy(ny) = 0;

% Variables for creating the U, k and ε equations
as = zeros(ny, 1);
ap = zeros(ny, 1);
an = zeros(ny, 1);
su = zeros(ny, 1);
sp = zeros(ny, 1);

iteration = 0; % Iteration count
while (max(res) > residue_limit)
    iteration = iteration + 1;
    % Solving for eddy-viscosity using the prandtl mixing length for some
    % initial iterations (say 100)
    if (iteration < 100)
        for i=2:ny-1
            nu_t(i) = (kappa * y_dns(i))^2 * (U(i+1)-U(i-1))/(dy(i-1)+dy(i));
        end
        nu_t(ny) = nu_t(ny-1);
    % Solving for eddy-viscosity using k-ε model
    else
        for i=2:ny-1
            nu_t(i) = c_mu * k(i)^2 / eps(i);
        end
        nu_t(ny) = nu_t(ny-1);
    end
    
    % Calculation of νₜ (Eddy viscosity) on the cell faces
    nu_t_surface = interpolation(nu_t, y_surface, y_dns, 1, ny-1);

    U_old = U;
    % Creating the U equations for the grid
    for i=2:ny-1
        if (i == 2)
            an(i) = (nu + nu_t_surface(i)) / dy(i);
            sp(i) = -(nu + nu_t_surface(i-1)) / dy(i-1);
            ap(i) = an(i) - sp(i);
            su(i) = (y_surface(i)-y_surface(i-1)) + (nu + nu_t_surface(i-1)) * U(1) / dy(i-1);
        elseif (i == ny-1)
            as(i) = (nu + nu_t_surface(i-1)) / dy(i-1);
            ap(i) = as(i);
            su(i) = y_surface(i)-y_surface(i-1);
        else
            as(i) = (nu + nu_t_surface(i-1)) / dy(i-1);
            an(i) = (nu + nu_t_surface(i)) / dy(i);
            ap(i) = as(i) + an(i);
            su(i) = y_surface(i)-y_surface(i-1);
        end
    end

    % Calculating the U
    U = gauss_seidel(as, ap, an, su, U, 2, ny-1);
    U(ny) = U(ny-1);

    % Calculating U using relaxation parameter
    U = (1-urf)*U_old + urf*U;

    % Calculating dUdy (Velocity Gradient)
    for i=1:ny-1
        if (i == 1)
            dUdy(i) = (U(i+1) - U(i)) / dy(i);
        else
            dUdy(i) = (U(i+1) - U(i-1)) / (dy(i-1) + dy(i));
        end
    end

    % Calculating friction velocity
    u_star = (nu*dUdy(1))^0.5;

    % Calculating yᐩ
    for i=1:ny
        y_plus(i) = u_star*y_dns(i)/nu;
    end

    % Calculating Pk (Production Rate)
    for i=1:ny-1
        Pk(i) = nu_t(i) * dUdy(i)^2;
    end

    % Creating the k equations for the grid
    k_old = k;
    for i=2:ny-1
        if (i == 2)
            an(i) = (nu + nu_t_surface(i)/sigma_k) / dy(i);
            sp(i) = -(nu + nu_t_surface(i-1)/sigma_k) / dy(i-1) - eps(i)*(y_surface(i)-y_surface(i-1))/k(i);
            ap(i) = an(i) - sp(i);
            su(i) = (y_surface(i)-y_surface(i-1))*Pk(i) + (nu + nu_t_surface(i-1)/sigma_k) * k(1) / dy(i-1);
        elseif (i == ny-1)
            as(i) = (nu + nu_t_surface(i-1)/sigma_k) / dy(i-1);
            sp(i) = -eps(i)*(y_surface(i)-y_surface(i-1))/k(i);
            ap(i) = as(i) - sp(i);
            su(i) = (y_surface(i)-y_surface(i-1))*Pk(i);
        else
            as(i) = (nu + nu_t_surface(i-1)/sigma_k) / dy(i-1);
            sp(i) = -eps(i)*(y_surface(i)-y_surface(i-1))/k(i);
            an(i) = (nu + nu_t_surface(i)/sigma_k) / dy(i);
            ap(i) = as(i) + an(i) - sp(i);
            su(i) = (y_surface(i)-y_surface(i-1))*Pk(i);
        end
    end

    % Calculating the k
    k = gauss_seidel(as, ap, an, su, k, 2, ny-1);
    k(ny) = k(ny-1);

    % Calculating k using relaxation parameter
    k = (1-urf)*k_old + urf*k;

    % Creating the ε equations for the grid
    eps_old = eps;
    for i=2:ny-1
        if (i == 2)
            an(i) = (nu + nu_t_surface(i)/sigma_eps) / dy(i);
            sp(i) = -c2_eps*eps(i)*(y_surface(i)-y_surface(i-1))/k(i);
            ap(i) = an(i) - sp(i);
            su(i) = c1_eps*Pk(i)*eps(i)*(y_surface(i)-y_surface(i-1))/k(i);
        elseif (i == ny-1)
            as(i) = (nu + nu_t_surface(i-1)/sigma_eps) / dy(i-1);
            sp(i) = -c2_eps*eps(i)*(y_surface(i)-y_surface(i-1))/k(i);
            ap(i) = as(i) - sp(i);
            su(i) = c1_eps*Pk(i)*eps(i)*(y_surface(i)-y_surface(i-1))/k(i);
        else
            as(i) = (nu + nu_t_surface(i-1)/sigma_eps) / dy(i-1);
            an(i) = (nu + nu_t_surface(i)/sigma_eps) / dy(i);
            sp(i) = -c2_eps*eps(i)*(y_surface(i)-y_surface(i-1)) / k(i);
            ap(i) = as(i) + an(i) - sp(i);
            su(i) = c1_eps*Pk(i)*eps(i)*(y_surface(i)-y_surface(i-1)) / k(i);
        end
    end

    % Checking for yᐩ ≤ 3 
    m = 0;
    for i=1:ny
        if (y_plus(i) > 3)
            m = i-1;
            break
        end
    end

    % Enforcing the given ε boundary condition at yᐩ ≤ 3
    for i=2:m
        eps(i) = 2 * nu * k(i) / y_dns(i)^2;
    end

    % Calculating the ε
    eps = gauss_seidel(as, ap, an, su, eps, m+1, ny-1);
    eps(1) = eps(2);
    eps(ny) = eps(ny-1);

    % Calculating ε using relaxation parameter
    eps = (1-urf)*eps_old + urf*eps;

    % Calculating the residuals
    res(1) = norm(U_old - U);
    res(2) = norm(k_old - k);
    res(3) = norm(eps_old - eps);

    % Printing the residuals
    fprintf(file, 'Iteration: %d\n', iteration);
    fprintf(file, 'U residual: %.4f\n', res(1));
    fprintf(file, 'k residual: %.4f\n', res(2));
    fprintf(file, 'ε residual: %.4f\n\n', res(3));
    fprintf('Iteration: %d\n', iteration);
    fprintf('U residual: %.4f\n', res(1));
    fprintf('k residual: %.4f\n', res(2));
    fprintf('ε residual: %.4f\n\n', res(3));
end

fprintf(file, 'The solution has converged!\n');
fprintf(file, 'τ = %.4f\n', nu*dUdy(1));
fprintf('The solution has converged!\n');
fprintf('τ = %.4f\n', nu*dUdy(1));
fprintf('Plotting the results...\n');

% Calculating the turbulence shear stress
tss = zeros(ny, 1); % Turbulence shear stress
for i=1:ny
    tss(i) = nu_t(i) * dUdy(i);
end

% Calculating turbulence diffusion of k
tdk = zeros(ny, 1);
tdk1 = zeros(ny, 1); % Turbulence diffusion of k
for i=2:ny-1
    tdk(i) = nu_t(i)*(k(i+1)-k(i-1))/((dy(i-1)+dy(i))*sigma_k);
end
for i=2:ny-1
    tdk1(i) = (tdk(i+1)-tdk(i-1))/(dy(i-1)+dy(i));
end

% Calculating viscous diffusion of k
vdk = zeros(ny, 1);
vdk1 = zeros(ny, 1); % Viscous diffusion of k
for i=1:ny-1
    if (i == 1)
        vdk(i) = nu*(k(i+1)-k(i))/dy(i);
    else
        r = dy(i) / dy(i-1);
        vdk(i) = nu*(k(i+1)-k(i-1))/(dy(i-1)+dy(i));
    end
end
for i=1:ny-1
    if (i == 1)
        vdk1(i) = (vdk(i+1)-vdk(i))/dy(i);
    else
        vdk1(i) = (vdk(i+1)-vdk(i-1))/(dy(i-1)+dy(i));
    end
end

% Calculating the Turbulence kinetic energy (k) using the DNS data
k_dns=0.5*(u2_dns+v2_dns+w2_dns);
eps_dns=dns_data(:,2)*ustar^4/nu; % eps is normalized by ustar^4/nu

% Saving the results
for i=1:ny
    fprintf(file1, ['%.4f    %.4f    %.4f    %.4f    %.4f    %.4f    %.4f    %.4f    ' ...
        '%.4f    %.4f    %.4f    %.4f    %.4f    %.4f    %.4f\n'], y_dns(i), u_dns(i), ...
        U(i), k_dns(i), k(i), eps_dns(i), eps(i), -uv_dns(i), tss(i), dns_data(i,3)/nu, Pk(i), dns_data(i,5)/nu, ...
        tdk1(i), dns_data(i, 6)/nu, vdk1(i));
end

% Closing the files
fclose(file);
fclose(file1);

% Plots for various parameters
figure(1)
plot(u_dns, y_dns, '-bo', U, y_dns, '-rx');
xlabel('U');
ylabel('y/h');
title('U-velocity');
legend('DNS', 'k-ε', 'Location', 'northwest');
legend boxoff;
saveas(1, './Plots/k_epsilon/U_velocity_plot', image_type);

figure(2)
plot(y_dns, k_dns, '-bo', y_dns, k, '-rx');
xlabel('y/h');
ylabel('k');
title('Turbulence kinetic energy');
legend('DNS', 'k-ε', 'Location', 'northeast');
legend boxoff;
saveas(2, './Plots/k_epsilon/TKE_plot', image_type);

figure(3)
plot(y_dns, eps_dns, '-bo', y_dns, eps, '-rx')
xlabel('y/h');
ylabel('\epsilon');
title('Dissipation rate of k');
legend('DNS', 'k-ε', 'Location', 'northeast');
legend boxoff
saveas(3, './Plots/k_epsilon/Dissipation_rate_of_k_plot', image_type);

figure(4)
plot(y_dns, -uv_dns, '-bo', y_dns, tss, '-rx');
xlabel('y/h');
ylabel('-<uv>');
title('Turbulence shear stress');
legend('DNS', 'k-ε', 'Location', 'northeast');
legend boxoff
saveas(4, './Plots/k_epsilon/Turbulent_shear_stress_plot', image_type);

figure(5)
plot(y_dns, dns_data(:,3)/nu, '-bo', y_dns, Pk, '-rx');
xlabel('y/h');
ylabel('P_k');
title('Production rate of k');
legend('DNS', 'k-ε', 'Location', 'northeast');
legend boxoff
saveas(5, './Plots/k_epsilon/Pk_plot', image_type);

figure(6)
plot(y_dns, dns_data(:,5)/nu, '-bo', y_dns, tdk1, '-rx');
xlabel('y/h');
title('Turbulent diffusion of k');
legend('DNS', 'k-ε', 'Location', 'northeast');
legend boxoff
saveas(6, './Plots/k_epsilon/Turbulent_diffusion_of_k_plot', image_type);

figure(7)
plot(y_dns, dns_data(:, 6)/nu, '-bo', y_dns, vdk1, '-rx');
xlabel('y/h');
title('Viscous diffusion of k');
legend('DNS', 'k-ε', 'Location', 'northeast');
legend boxoff
saveas(7, './Plots/k_epsilon/Viscous_diffusion_of_k_plot', image_type);

figure(8)
plot(dns_data(:, 1), -eps_dns, '-*c', dns_data(:, 1), dns_data(:,3)/nu, '--k', ...
    dns_data(:, 1), dns_data(:,5)/nu, '-xb', dns_data(:, 1), dns_data(:, 6)/nu, '-or');
xlabel('yᐩ');
title('Budget for k (DNS)');
legend('Dissipation rate of k', 'Production rate of k', ...
    'Turbulent diffusion of k', 'Viscous diffusion of k', 'Location', 'northeast');
legend boxoff
saveas(8, './Plots/k_epsilon/Budget for Turbulence (DNS)', image_type);

figure(9)
plot(y_plus, -eps, '-*c', y_plus, Pk, '--k', ...
    y_plus, tdk1, '-xb', y_plus, vdk1, '-or');
xlabel('yᐩ');
title('Budget for k (k-ε model)');
legend('Dissipation rate of k', 'Production rate of k', ...
    'Turbulent diffusion of k', 'Viscous diffusion of k', 'Location', 'northeast');
legend boxoff
saveas(9, './Plots/k_epsilon/Budget for Turbulence (k-epsilon)', image_type);

% Function for solving system of equations using Gauss-Siedel method
function x = gauss_seidel(as, ap, an, su, x, x_initial, x_final)
    res_limit = 1e-4; % Convergence criteria
    x_res = 2*res_limit; % Initial residue 
    while (x_res > res_limit)
        x_old = x;
        for i=x_initial:x_final
            x(i) = ((su(i)+as(i)*x(i-1)+an(i)*x(i+1))/ap(i));
        end
        x_res = norm(x-x_old); % Norm calculation
    end
end

% Function for interpolation of νₜ (Eddy viscosity)
function nu_t_surface = interpolation(nu_t, x_surface, x_node, x_initial, x_final)
    nu_t_surface(x_initial) = nu_t(x_initial);
    nu_t_surface(x_final) = nu_t(x_final);
    for i=x_initial+1:x_final-1
        nu_t_surface(i) = nu_t(i)+(nu_t(i+1)-nu_t(i))*(x_surface(i)-x_node(i))/...
        (x_node(i+1)-x_node(i));
    end
end
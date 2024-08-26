% Matlab file for plotting the comparison plots

clc; close all; clear;
image_type = 'epsc';

% Simulated data
k_eps = load('simulated_data_k_epsilon.dat');
rng_k_eps = load('simulated_data_rng_k_epsilon.dat');
k_omega = load('simulated_data_k_omega.dat');

% Plots for various parameters
figure(1)
plot(k_eps(:,2), k_eps(:,1), '-bo', k_eps(:,3), k_eps(:,1), '-g+', rng_k_eps(:,3), k_eps(:,1), '-r*' ...
    , k_omega(:,3), k_eps(:,1), '-kv');
xlabel('U');
ylabel('y/h');
title('U-velocity');
legend('DNS', 'k-ε', 'RNG k-ε', 'k-ω', 'Location', 'northwest');
legend boxoff;
saveas(1, './Plots/comparison/U_velocity_comparison_plot', image_type);

figure(2)
plot(k_eps(:,1), k_eps(:,4), '-bo', k_eps(:,1), k_eps(:,5), '-g+', k_eps(:,1), rng_k_eps(:,5), '-r*' ...
    , k_eps(:,1), k_omega(:,5), '-kv');
xlabel('y/h');
ylabel('k');
title('Turbulence kinetic energy');
legend('DNS', 'k-ε', 'RNG k-ε', 'k-ω', 'Location', 'northeast');
legend boxoff;
saveas(2, './Plots/comparison/TKE_comparison_plot', image_type);

figure(3)
plot(k_eps(:,1), k_eps(:,6), '-bo', k_eps(:,1), k_eps(:,7), '-g+', k_eps(:,1), rng_k_eps(:,7), '-r*' ...
    , k_eps(:,1), k_omega(:,7), '-kv');
xlabel('y/h');
ylabel('\epsilon');
title('Dissipation rate of k');
legend('DNS', 'k-ε', 'RNG k-ε', 'k-ω', 'Location', 'northeast');
legend boxoff
saveas(3, './Plots/comparison/Dissipation_rate_of_k_comparison_plot', image_type);

figure(4)
plot(k_eps(:,1), k_eps(:,8), '-bo', k_eps(:,1), k_eps(:,9), '-g+', k_eps(:,1), rng_k_eps(:,9), '-r*' ...
    , k_eps(:,1), k_omega(:,9), '-kv');
xlabel('y/h');
ylabel('-<uv>');
title('Turbulence shear stress');
legend('DNS', 'k-ε', 'RNG k-ε', 'k-ω', 'Location', 'northeast');
legend boxoff
saveas(4, './Plots/comparison/Turbulent_shear_stress_comparison_plot', image_type);

figure(5)
plot(k_eps(:,1), k_eps(:,10), '-bo', k_eps(:,1), k_eps(:,11), '-g+', k_eps(:,1), rng_k_eps(:,11), '-r*' ...
    , k_eps(:,1), k_omega(:,11), '-kv');
xlabel('y/h');
ylabel('P_k');
title('Production rate of k');
legend('DNS', 'k-ε', 'RNG k-ε', 'k-ω', 'Location', 'northeast');
legend boxoff
saveas(5, './Plots/comparison/Pk_comparison_plot', image_type);

figure(6)
plot(k_eps(:,1), k_eps(:,12), '-bo', k_eps(:,1), k_eps(:,13), '-g+', k_eps(:,1), rng_k_eps(:,13), '-r*' ...
    , k_eps(:,1), k_omega(:,13), '-kv');
xlabel('y/h');
title('Turbulent diffusion of k');
legend('DNS', 'k-ε', 'RNG k-ε', 'k-ω', 'Location', 'northeast');
legend boxoff
saveas(6, './Plots/comparison/Turbulent_diffusion_of_k_comparison_plot', image_type);

figure(7)
plot(k_eps(:,1), k_eps(:,14), '-bo', k_eps(:,1), k_eps(:,15), '-g+', k_eps(:,1), rng_k_eps(:,15), '-r*' ...
    , k_eps(:,1), k_omega(:,15), '-kv');
xlabel('y/h');
title('Viscous diffusion of k');
legend('DNS', 'k-ε', 'RNG k-ε', 'k-ω', 'Location', 'northeast');
legend boxoff
saveas(7, './Plots/comparison/Viscous_diffusion_of_k_comparison_plot', image_type);
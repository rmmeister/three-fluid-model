%% INITIALIZE MATLAB
clear all
clc
close all
format long
%% Solving Y with Levenberg-Marquardt Optimization
Y0 = [0.949; 0.05; .001; 30; 2; 25];     % initial guess for unknowns vector
alphaGvG = 10:1:50;                      % product of gas velocity and void fraction
Y = zeros(6, length(alphaGvG));          % unknowns vector

fprintf('Solving for vector of unknowns at %g points...\n', length(alphaGvG));
for i = 1:length(alphaGvG)
    [Y(:, i), residualVal, yFlag, output] = fsolve(@(Y) calcResiduals(Y, alphaGvG(i)), Y0, ...
        optimoptions('fsolve','Display','off'));  
    [ F, GammaDL(i), GammaLD(i), dD(i), RGD(i), RLW(i), RGL(i), RDL(i) ] = calcResiduals(  Y(:, i), alphaGvG(i) );
end
if yFlag >= 1
    fprintf('\nSimulation converged successfully.\n');
end

%% RESULTS
gasDensity = 1.725; % Kg/m3
liquidDensity = 998; % Kg/m3

% PLOTING THE RESULTS
col = [0.850980401039124 0.325490206480026 0.0980392172932625];
figure('Position', [10 10 1400 700]) 
subplot(231);
semilogy(alphaGvG, Y(1, :), 'Color', 'b' ,'LineWidth', 3);
grid on; grid minor;
xlabel('\alpha_G v_G, m/s');
ylabel('\alpha');
title('Gas, liquid film and droplet fractions');
hold on 
semilogy(alphaGvG, Y(2, :), 'Color', 'r', 'LineWidth', 3); 
hold on
semilogy(alphaGvG, Y(3, :), 'Color', 'g', 'LineWidth', 3);
legend('\alpha_G', '\alpha_L', '\alpha_D', 'Location', 'NorthEast');
xlim([0 60]);
hold off

subplot(232);
plot(alphaGvG, gasDensity.*Y(1, :).*Y(4, :), 'Color', 'b', 'LineWidth', 3); 
xlabel('\alpha_G v_G, m/s');
ylabel('\alpha v \rho, kg/m^2/s');
title('Gas, liquid film and droplet superficial mass flow rates');
hold on
plot(alphaGvG, liquidDensity.*Y(2, :).*Y(5, :), 'Color', 'r', 'LineWidth', 3);
hold on
plot(alphaGvG, liquidDensity.*Y(3, :).*Y(6, :), 'Color', 'g', 'LineWidth', 3); grid on;
legend('\alpha_G v_G \rho_G', '\alpha_L v_L \rho_L', '\alpha_D v_D \rho_L', 'Location', 'NorthEast');
xlim([0 60]);
ylim([0 100]);
hold off

subplot(233);
plot(alphaGvG, Y(4, :), 'Color', 'b', 'LineWidth', 3); grid on;
xlabel('\alpha_G v_G, m/s');
ylabel('v, m/s');
title('Gas, liquid film and droplet velocities');
hold on
plot(alphaGvG, Y(5, :), 'Color', 'r' ,'LineWidth', 3);
hold on 
plot(alphaGvG, Y(6, :), 'Color', 'g', 'LineWidth', 3); 
legend('v_G', 'v_L', 'v_D', 'Location', 'NorthWest');
xlim([0 60]);
hold off


subplot(234);
plot(alphaGvG, GammaDL, 'Color', col, 'LineWidth', 3); grid on;
xlabel('\alpha_G v_G, m/s');
ylabel('\Gamma_{LDi}, \Gamma_{DLi}, kg/m^3/s');
title('Liquid film entrainment and droplet deposition');
xlim([0 60]);
hold on 
plot(alphaGvG, real(GammaLD), 'Color', 'b', 'LineWidth', 3);
ylim([0 20]);
hold off

subplot(235);
plot(alphaGvG, dD, 'Color', col, 'LineWidth', 3); grid on; 
xlabel('\alpha_G v_G, m/s');
ylabel('d_D, m');
title('Mean droplet (Sauter) diameter');
xlim([0 60]);

subplot(236);
semilogy(alphaGvG, RGL, 'Color', 'g', 'LineWidth', 3); grid on; grid minor; hold all
semilogy(alphaGvG, RGD, 'Color', 'm', 'LineWidth', 3); 
semilogy(alphaGvG, abs(RLW), 'Color', 'y', 'LineWidth', 3); 
semilogy(alphaGvG, RDL, 'Color', 'c', 'LineWidth', 3); 
semilogy(alphaGvG, abs(GammaDL.*Y(6, :)), 'Color', col, 'LineWidth', 3); 
semilogy(alphaGvG, real(GammaLD).*abs(Y(5, :)), 'Color', 'b', 'LineWidth', 3); hold off
title('Volume-specific Forces, kg/m^2/s^2');
legend('R_{GL}', 'R_{GD}', 'R_{LW}', 'R_{DL}', '\Gamma_{DL} v_D', '\Gamma_{LD} v_L', 'Location', ...
    'SouthEast');
xlabel('\alpha_G v_G, m/s');
xlim([0 60]);

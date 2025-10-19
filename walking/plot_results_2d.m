function plot_results_2d(results)
% PLOT_RESULTS_2D - Visualization of simulation results
%
% Creates comprehensive plots of:
%   - Satellite attitude (critical for disturbance analysis)
%   - Momentum conservation
%   - Reaction wheel usage
%   - Walker trajectory
%   - Contact forces
%
% Input:
%   results - Results structure from simulate_experiment_2d

fprintf('Generating plots...\n');

%% Create figure
fig = figure('Position', [100, 100, 1400, 900], 'Name', '2D Orbital Walker Results');

%% Subplot 1: Satellite Attitude (CRITICAL)
subplot(3, 3, 1);
plot(results.t, rad2deg(results.satellite.theta), 'b-', 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Satellite Rotation [deg]');
title('Satellite Attitude θ_{sat}');
yline(0, 'k--', 'LineWidth', 1);

% Highlight max deviation
[max_theta, idx_max] = max(abs(results.satellite.theta));
hold on;
plot(results.t(idx_max), rad2deg(results.satellite.theta(idx_max)), 'ro', ...
     'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(results.t(idx_max), rad2deg(results.satellite.theta(idx_max)), ...
     sprintf('  Max: %.4f°', rad2deg(max_theta)), 'FontSize', 9);

%% Subplot 2: Momentum Conservation
subplot(3, 3, 2);
plot(results.t, results.momentum.conservation_error, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('||h_{total} - h_{init}|| [N·m·s]');
title('Momentum Conservation Error');
set(gca, 'YScale', 'log');

%% Subplot 3: Reaction Wheel Usage
subplot(3, 3, 3);
h_max = 30;  % From environment setup
plot(results.t, results.RW.h, 'b-', 'LineWidth', 2);
hold on;
yline(h_max, 'r--', 'LineWidth', 1.5, 'h_{max}');
yline(-h_max, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('RW Momentum [N·m·s]');
title('Reaction Wheel Storage');
ylim([-h_max*1.2, h_max*1.2]);

%% Subplot 4: Walker Trajectory (Top View)
subplot(3, 3, 4);
plot(results.walker.pos(1,:), results.walker.pos(3,:), 'b-', 'LineWidth', 2);
hold on;
plot(results.walker.pos(1,1), results.walker.pos(3,1), 'go', ...
     'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(results.walker.pos(1,end), results.walker.pos(3,end), 'ro', ...
     'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'End');

% Plot satellite
plot(results.satellite.pos(1,:), results.satellite.pos(3,:), 'k--', 'LineWidth', 1);
grid on;
xlabel('X [m]');
ylabel('Z [m]');
title('Trajectory (Top View)');
legend('Walker', 'Start', 'End', 'Satellite');
axis equal;

%% Subplot 5: Walker Torso Angle
subplot(3, 3, 5);
plot(results.t, rad2deg(results.walker.theta), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Walker Rotation [deg]');
title('Walker Torso Orientation θ_{walker}');

%% Subplot 6: Joint Angles
subplot(3, 3, 6);
hold on;
for i = 1:6
    plot(results.t, rad2deg(results.walker.q_joints(i,:)), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Joint %d', i));
end
grid on;
xlabel('Time [s]');
ylabel('Joint Angles [deg]');
title('Walker Joint Angles');
legend('Location', 'best');

%% Subplot 7: Contact Forces (if available)
subplot(3, 3, 7);
if isfield(results.contacts, 'lambda') && any(results.contacts.lambda(:))
    lambda_norms = vecnorm(results.contacts.lambda);
    plot(results.t, lambda_norms, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('||λ|| [N]');
    title('Contact Force Magnitude');
else
    text(0.5, 0.5, 'No contact force data', ...
         'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis off;
end

%% Subplot 8: Constraint Violation
subplot(3, 3, 8);
if isfield(results.contacts, 'violation')
    semilogy(results.t, results.contacts.violation, 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('||Φ||');
    title('Constraint Violation');
else
    text(0.5, 0.5, 'No constraint data', ...
         'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis off;
end

%% Subplot 9: Phase Portrait (Satellite)
subplot(3, 3, 9);
% Compute angular velocity from quaternion derivative (simplified)
omega_sat = gradient(results.satellite.theta) / mean(diff(results.t));
plot(rad2deg(results.satellite.theta), rad2deg(omega_sat), 'b-', 'LineWidth', 1.5);
hold on;
plot(rad2deg(results.satellite.theta(1)), rad2deg(omega_sat(1)), 'go', ...
     'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(rad2deg(results.satellite.theta(end)), rad2deg(omega_sat(end)), 'ro', ...
     'MarkerSize', 10, 'MarkerFaceColor', 'r');
grid on;
xlabel('θ_{sat} [deg]');
ylabel('ω_{sat} [deg/s]');
title('Satellite Phase Portrait');

%% Add overall title
sgtitle('2D Orbital Walker Simulation Results', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('✓ Plots generated\n');

%% Print key metrics
fprintf('\n=== Key Performance Metrics ===\n');
fprintf('Satellite attitude deviation: %.4f deg (%.2e rad)\n', ...
        rad2deg(results.stats.max_satellite_rotation), results.stats.max_satellite_rotation);
fprintf('RW utilization: %.1f%% (%.2f / %.0f N·m·s)\n', ...
        100 * results.stats.max_RW_momentum / 30, results.stats.max_RW_momentum, 30);
fprintf('Momentum conservation: %.2e N·m·s\n', results.stats.max_conservation_error);

if results.stats.max_satellite_rotation < deg2rad(0.1)
    fprintf('✓ EXCELLENT: Satellite disturbance < 0.1 deg\n');
elseif results.stats.max_satellite_rotation < deg2rad(1)
    fprintf('✓ GOOD: Satellite disturbance < 1 deg\n');
else
    fprintf('⚠ WARNING: Satellite disturbance > 1 deg\n');
end

end
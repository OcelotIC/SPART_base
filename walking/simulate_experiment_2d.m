function results = simulate_experiment_2d(traj_type)
% SIMULATE_EXPERIMENT_2D - Main simulation loop for 2D orbital walker
%
% Inputs:
%   traj_type - Trajectory type: 'line' (default), 'circle', 'static'
%
% Outputs:
%   results - Structure with simulation results

if nargin < 1
    traj_type = 'line';
end

fprintf('=== 2D Orbital Walker Simulator ===\n\n');

%% 1. SETUP
fprintf('[1/6] Setting up models...\n');

% Create robot models
[robot_walker, ~] = create_walker_2d();
[robot_sat, satellite_params] = create_satellite_2d();

% Environment parameters
env = setup_environment_2d();
env.satellite_params = satellite_params;

%% 2. GENERATE TRAJECTORY
fprintf('[2/6] Generating trajectory...\n');
traj = generate_trajectory_2d(env, traj_type);
N = length(traj.t);

%% 3. INITIALIZE STATE
fprintf('[3/6] Initializing state...\n');
state = initialize_state_2d(env, robot_walker, robot_sat, satellite_params);

%% 4. PRE-ALLOCATE RESULTS
fprintf('[4/6] Allocating memory...\n');
results = struct();
results.t = traj.t;
results.walker = struct();
results.satellite = struct();
results.RW = struct();
results.momentum = struct();
results.contacts = struct();

% Pre-allocate arrays
results.walker.pos = zeros(3, N);
results.walker.theta = zeros(1, N);
results.walker.q_joints = zeros(6, N);
results.satellite.pos = zeros(3, N);
results.satellite.theta = zeros(1, N);
results.RW.h = zeros(1, N);
results.momentum.h_total = zeros(3, N);
results.momentum.conservation_error = zeros(1, N);
results.contacts.lambda = zeros(6, N);
results.contacts.violation = zeros(1, N);

%% 5. SIMULATION LOOP
fprintf('[5/6] Running simulation...\n');
tic;

for k = 1:N
    % Progress indicator
    if mod(k, floor(N/10)) == 0
        fprintf('  Progress: %d%%\n', round(100*k/N));
    end
    
    % Compute system dynamics
    sys = compute_coupled_dynamics_2d(state, robot_walker, robot_sat, satellite_params, env);
    
    % Controller (PD on joints)
    q_desired = traj.q_joints(:, k);
    qd_desired = zeros(6, 1);  % Zero velocity for static pose
    tau_joints = pd_controller_2d(state, q_desired, qd_desired, env);
    
    % Integrate dynamics
    if k < N
        state = integrate_dynamics_2d(state, tau_joints, sys, env);
    end
    
    % Save results
    results.walker.pos(:, k) = state.r0_walker;
    R_walker = quat_DCM(state.quat_walker);
    results.walker.theta(k) = atan2(R_walker(3,1), R_walker(1,1));
    results.walker.q_joints(:, k) = state.qm_walker;
    
    results.satellite.pos(:, k) = state.r0_sat;
    R_sat = quat_DCM(state.quat_sat);
    results.satellite.theta(k) = atan2(R_sat(3,1), R_sat(1,1));
    
    results.RW.h(k) = state.h_RW;
    results.momentum.h_total(:, k) = sys.h_total;
    results.momentum.conservation_error(k) = norm(sys.h_total - state.h_total_initial);
    
    if isfield(state, 'lambda')
        results.contacts.lambda(:, k) = state.lambda;
        results.contacts.violation(k) = state.constraint_violation;
    end
end

sim_time = toc;
fprintf('  Simulation completed in %.2f s (%.1fx real-time)\n', ...
        sim_time, traj.t(end)/sim_time);

%% 6. POST-PROCESSING
fprintf('[6/6] Post-processing...\n');

% Compute statistics
results.stats = struct();
results.stats.max_satellite_rotation = max(abs(results.satellite.theta));
results.stats.max_RW_momentum = max(abs(results.RW.h));
results.stats.max_conservation_error = max(results.momentum.conservation_error);
results.stats.max_constraint_violation = max(results.contacts.violation);

fprintf('\n=== Simulation Summary ===\n');
fprintf('Duration: %.1f s (%d steps)\n', traj.t(end), N);
fprintf('Max satellite rotation: %.4f deg\n', rad2deg(results.stats.max_satellite_rotation));
fprintf('Max RW momentum: %.2f / %.0f N·m·s\n', ...
        results.stats.max_RW_momentum, satellite_params.RW.h_max);
fprintf('Max momentum conservation error: %.2e N·m·s\n', results.stats.max_conservation_error);
fprintf('Max constraint violation: %.2e\n', results.stats.max_constraint_violation);

fprintf('\n✓ Simulation complete!\n\n');

end

%% Helper Function: Initialize State
function state = initialize_state_2d(env, robot_walker, robot_sat, satellite_params)
    % Initialize complete system state
    
    state = struct();
    
    % Walker
    state.r0_walker = env.initial.walker.position;
    theta_w = env.initial.walker.orientation;
    state.quat_walker = [0; sin(theta_w/2); 0; cos(theta_w/2)];  % Rotation around Y
    state.qm_walker = env.initial.walker.joints;
    state.u0_walker = zeros(6, 1);  % [ωx; ωy; ωz; vx; vy; vz]
    state.um_walker = zeros(6, 1);
    
    % Satellite
    state.r0_sat = env.initial.satellite.position;
    theta_s = env.initial.satellite.orientation;
    state.quat_sat = [0; sin(theta_s/2); 0; cos(theta_s/2)];
    state.u0_sat = zeros(6, 1);
    
    % Reaction Wheel
    state.h_RW = env.initial.RW.h_stored;
    
    % Initial momentum (for conservation check)
    state.h_total_initial = [0; 0; 0];
    
    fprintf('  Initial state set\n');
end
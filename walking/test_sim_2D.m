%% TEST_SIMULATOR_2D - Main test script for 2D orbital walker simulator
%
% This script demonstrates the complete simulation workflow:
%   1. Model creation
%   2. Trajectory generation
%   3. Dynamics simulation
%   4. Results visualization
%
% Usage:
%   Simply run this script: >> test_simulator_2d
%
% Requirements:
%   - SPART toolbox must be in MATLAB path
%   - All simulator_2d/ functions must be available

clear; clc; close all;

fprintf('\n');
fprintf('╔════════════════════════════════════════════╗\n');
fprintf('║   2D ORBITAL WALKER SIMULATOR TEST        ║\n');
fprintf('║   Testing constrained dynamics solver     ║\n');
fprintf('╚════════════════════════════════════════════╝\n');
fprintf('\n');

%% SETUP
fprintf('┌─ Setup Phase ────────────────────────────┐\n');

% Add SPART to path (adjust if needed)
if exist('Kinematics.m', 'file') ~= 2
    fprintf('│ ⚠ SPART not found in path              │\n');
    fprintf('│   Add SPART using: addpath(''path/to/SPART/src'') │\n');
    fprintf('└──────────────────────────────────────────┘\n');
    return;
else
    fprintf('│ ✓ SPART toolbox found                   │\n');
end

fprintf('└──────────────────────────────────────────┘\n\n');

%% TEST 1: Model Creation
fprintf('┌─ Test 1: Model Creation ─────────────────┐\n');
try
    [robot_walker, ~] = create_walker_2d();
    fprintf('│ ✓ Walker model created                  │\n');
    
    [robot_sat, sat_params] = create_satellite_2d();
    fprintf('│ ✓ Satellite model created               │\n');
    
    env = setup_environment_2d();
    fprintf('│ ✓ Environment configured                │\n');
    
    fprintf('└──────────────────────────────────────────┘\n\n');
catch ME
    fprintf('│ ✗ ERROR: %s\n', ME.message);
    fprintf('└──────────────────────────────────────────┘\n');
    return;
end

%% TEST 2: Trajectory Generation
fprintf('┌─ Test 2: Trajectory Generation ──────────┐\n');
try
    traj = generate_trajectory_2d(env, 'line');
    fprintf('│ ✓ Trajectory generated (%d points)      │\n', length(traj.t));
    
    fprintf('└──────────────────────────────────────────┘\n\n');
catch ME
    fprintf('│ ✗ ERROR: %s\n', ME.message);
    fprintf('└──────────────────────────────────────────┘\n');
    return;
end

%% TEST 3: Single Step Dynamics
fprintf('┌─ Test 3: Single Step Dynamics ───────────┐\n');
try
    % Initialize state
    state = struct();
    state.r0_walker = [0; 0; 0.8];
    state.quat_walker = [0; 0; 0; 1];
    state.qm_walker = deg2rad([30; -60; 30; -30; -60; -30]);
    state.u0_walker = zeros(6, 1);
    state.um_walker = zeros(6, 1);
    state.r0_sat = [0; 0; 0];
    state.quat_sat = [0; 0; 0; 1];
    state.u0_sat = zeros(6, 1);
    state.h_RW = 0;
    state.h_total_initial = [0; 0; 0];
    
    % Add satellite params to env
    env.satellite_params = sat_params;
    
    % Compute dynamics
    sys = compute_coupled_dynamics_2d(state, robot_walker, robot_sat, sat_params, env);
    fprintf('│ ✓ Dynamics computed                     │\n');
    fprintf('│   M matrix: %dx%d                       │\n', size(sys.M_full));
    fprintf('│   Constraints: %d equations              │\n', length(sys.Phi));
    
    % Test integration
    tau_joints = zeros(6, 1);
    state_next = integrate_dynamics_2d(state, tau_joints, sys, env);
    fprintf('│ ✓ Time integration successful           │\n');
    
    fprintf('└──────────────────────────────────────────┘\n\n');
catch ME
    fprintf('│ ✗ ERROR: %s\n', ME.message);
    fprintf('│   %s (line %d)\n', ME.stack(1).file, ME.stack(1).line);
    fprintf('└──────────────────────────────────────────┘\n');
    return;
end

%% TEST 4: Full Simulation
fprintf('┌─ Test 4: Full Simulation ────────────────┐\n');
fprintf('│ Running simulation...                    │\n');
fprintf('│ (This may take 10-30 seconds)            │\n');
fprintf('└──────────────────────────────────────────┘\n\n');

try
    results = simulate_experiment_2d('line');
    fprintf('\n┌─ Simulation Results ──────────────────────┐\n');
    fprintf('│ ✓ Simulation completed successfully     │\n');
    fprintf('└──────────────────────────────────────────┘\n\n');
catch ME
    fprintf('\n┌─ Simulation Error ────────────────────────┐\n');
    fprintf('│ ✗ ERROR: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('│   File: %s\n', ME.stack(1).file);
        fprintf('│   Line: %d\n', ME.stack(1).line);
    end
    fprintf('└──────────────────────────────────────────┘\n');
    return;
end

%% TEST 5: Visualization
fprintf('┌─ Test 5: Visualization ───────────────────┐\n');
try
    plot_results_2d(results);
    fprintf('│ ✓ Plots generated                       │\n');
    fprintf('└──────────────────────────────────────────┘\n\n');
catch ME
    fprintf('│ ✗ ERROR: %s\n', ME.message);
    fprintf('└──────────────────────────────────────────┘\n');
    return;
end

%% FINAL SUMMARY
fprintf('╔════════════════════════════════════════════╗\n');
fprintf('║         ALL TESTS PASSED ✓                 ║\n');
fprintf('╚════════════════════════════════════════════╝\n\n');

fprintf('Next steps:\n');
fprintf('  1. Verify constraint satisfaction (should be < 1e-6)\n');
fprintf('  2. Check momentum conservation (should be < 1e-4)\n');
fprintf('  3. Adjust controller gains if needed\n');
fprintf('  4. Test different trajectories: ''circle'', ''static''\n');
fprintf('  5. Implement NMPC and QP controllers\n\n');

fprintf('To run again: >> test_simulator_2d\n');
fprintf('To run with different trajectory: >> results = simulate_experiment_2d(''circle'')\n\n');
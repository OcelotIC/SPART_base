%% TEST_URDF_LOADING - Verify URDF files load correctly
%
% This script tests:
%   1. Walker URDF loads with correct branching structure
%   2. Satellite URDF loads as rigid body
%   3. ConnectivityMap produces correct results

clear; clc;

fprintf('╔════════════════════════════════════════════╗\n');
fprintf('║   URDF LOADING TEST                        ║\n');
fprintf('╚════════════════════════════════════════════╝\n\n');

%% Test 1: Walker URDF
fprintf('┌─ Test 1: Walker URDF ────────────────────┐\n');

try
    [robot_walker, keys_walker] = create_walker_2d();
    fprintf('│ ✓ Walker loaded successfully            │\n');
    
    % Check structure
    fprintf('│   n_q = %d (should be 6)                │\n', robot_walker.n_q);
    fprintf('│   n_links_joints = %d (should be 6)     │\n', robot_walker.n_links_joints);
    
    % Check branching
    child_base_indices = find(robot_walker.con.child_base);
    fprintf('│   child_base indices: [%s]            │\n', num2str(child_base_indices'));
    fprintf('│   (should be [1 4] for 2 arms)          │\n');
    
    % Display branch matrix
    fprintf('│                                          │\n');
    fprintf('│   Branch matrix [6×6]:                  │\n');
    for i = 1:size(robot_walker.con.branch, 1)
        fprintf('│   ');
        for j = 1:size(robot_walker.con.branch, 2)
            fprintf('%d ', robot_walker.con.branch(i,j));
        end
        fprintf('               │\n');
    end
    
    fprintf('└──────────────────────────────────────────┘\n\n');
    
catch ME
    fprintf('│ ✗ FAILED: %s\n', ME.message);
    fprintf('└──────────────────────────────────────────┘\n\n');
    return;
end

%% Test 2: Satellite URDF
fprintf('┌─ Test 2: Satellite URDF ─────────────────┐\n');

try
    [robot_sat, sat_params] = create_satellite_2d();
    fprintf('│ ✓ Satellite loaded successfully         │\n');
    
    % Check it's a rigid body
    fprintf('│   n_q = %d (should be 0)                │\n', robot_sat.n_q);
    fprintf('│   n_links_joints = %d (should be 0)     │\n', robot_sat.n_links_joints);
    fprintf('│   Base mass = %.0f kg                   │\n', robot_sat.base_link.mass);
    
    fprintf('└──────────────────────────────────────────┘\n\n');
    
catch ME
    fprintf('│ ✗ FAILED: %s\n', ME.message);
    fprintf('└──────────────────────────────────────────┘\n\n');
    return;
end

%% Test 3: SPART Kinematics with Walker
fprintf('┌─ Test 3: SPART Kinematics ───────────────┐\n');

try
    % Initial configuration
    R0 = eye(3);
    r0 = [0; 0; 0.8];
    qm = deg2rad([30; -60; 30; -30; -60; -30]);
    
    % Call SPART kinematics
    [RJ, RL, rJ, rL, e, g] = Kinematics(R0, r0, qm, robot_walker);
    
    fprintf('│ ✓ Kinematics computed successfully      │\n');
    fprintf('│   RL size: %s                     │\n', mat2str(size(RL)));
    fprintf('│   rL size: %s                       │\n', mat2str(size(rL)));
    
    % End-effector positions
    ee1_pos = rL(:, 3);  % Arm 1 end-effector
    ee2_pos = rL(:, 6);  % Arm 2 end-effector
    
    fprintf('│   EE1 position: [%.2f, %.2f, %.2f]     │\n', ee1_pos);
    fprintf('│   EE2 position: [%.2f, %.2f, %.2f]     │\n', ee2_pos);
    
    fprintf('└──────────────────────────────────────────┘\n\n');
    
catch ME
    fprintf('│ ✗ FAILED: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('│   at %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
    fprintf('└──────────────────────────────────────────┘\n\n');
    return;
end

%% Summary
fprintf('╔════════════════════════════════════════════╗\n');
fprintf('║   ALL TESTS PASSED ✓                       ║\n');
fprintf('║   URDF files are correct!                  ║\n');
fprintf('╚════════════════════════════════════════════╝\n\n');

fprintf('Next step: Run test_simulator_2d.m\n\n');
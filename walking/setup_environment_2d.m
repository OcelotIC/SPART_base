function env = setup_environment_2d()
% SETUP_ENVIRONMENT_2D - Global parameters for 2D orbital walking simulation
%
% Output:
%   env - Structure containing all simulation parameters
%
% Author: Generated for orbital walking robot project

%% Simulation Parameters
env.sim = struct();
env.sim.dt = 0.01;           % Time step [s]
env.sim.T_final = 20.0;      % Total simulation time [s]
env.sim.gravity = 0;         % No gravity in orbit [m/s²]

%% Numerical Solver Parameters
env.solver = struct();
env.solver.method = 'implicit_euler';  % or 'rk4', 'ode45'
env.solver.constraint_tol = 1e-6;      % Contact constraint tolerance
env.solver.momentum_tol = 1e-4;        % Momentum conservation check
env.solver.max_iterations = 50;        % Max iterations for constraint solver

%% Contact Model Parameters
env.contact = struct();
env.contact.type = 'bilateral_rigid';  % Rigid docking (no compliance)
env.contact.n_contacts = 2;            % Double contact (both arms)

% Contact force limits (HOTDOCK-like interface)
env.contact.f_max = 500;     % Maximum normal/tangential force [N]
env.contact.tau_max = 50;    % Maximum contact torque [N·m]

% Stabilization (Baumgarte method)
env.contact.alpha_stab = 5.0;  % Damping coefficient
env.contact.beta_stab = 10.0;  % Stiffness coefficient

%% Momentum Conservation Constraint
env.momentum = struct();
env.momentum.h_total_initial = [0; 0; 0];  % Initial angular momentum [N·m·s]
env.momentum.enforce_conservation = true;  % Hard constraint vs soft monitoring

%% Physical Constants
env.physics = struct();
env.physics.frame = 'inertial';  % Reference frame
env.physics.dimension = '2D';    % XZ plane, rotation around Y

%% Initial Conditions (default, can be overridden)
env.initial = struct();

% Walker initial pose
env.initial.walker.position = [0; 0; 0.8];  % [x; y; z] torso position
env.initial.walker.orientation = 0;         % θ torso [rad]
env.initial.walker.joints = [
    deg2rad(30);   % arm1_joint1 (shoulder)
    deg2rad(-60);  % arm1_joint2 (elbow bent)
    deg2rad(30);   % arm1_joint3 (wrist)
    deg2rad(-30);  % arm2_joint1
    deg2rad(-60);  % arm2_joint2
    deg2rad(-30)   % arm2_joint3
];

% Satellite initial pose
env.initial.satellite.position = [0; 0; 0];  % [x; y; z] center
env.initial.satellite.orientation = 0;       % θ satellite [rad]
env.initial.satellite.velocity = [0; 0; 0];  % [vx; vy; vz]
env.initial.satellite.omega = [0; 0; 0];     % [ωx; ωy; ωz]

% Reaction wheel initial state
env.initial.RW.h_stored = 0;                 % Initial momentum [N·m·s]
env.initial.RW.omega = 0;                    % Initial wheel speed [rad/s]

%% Controller Parameters (for testing)
env.controller = struct();
env.controller.type = 'pd_joint';  % Simple PD controller
env.controller.Kp = 100 * eye(6);  % Proportional gain
env.controller.Kd = 20 * eye(6);   % Derivative gain

%% Data Logging
env.logging = struct();
env.logging.save_rate = 10;  % Save every N steps (to reduce memory)
env.logging.plot_realtime = false;
env.logging.verbose = true;

%% Visualization
env.viz = struct();
env.viz.animate = true;
env.viz.frame_rate = 30;  % Frames per second
env.viz.axis_limits = [-3, 3, -2, 3];  % [xmin, xmax, zmin, zmax]
env.viz.show_momentum = true;
env.viz.show_contacts = true;

%% Safety Limits
env.limits = struct();
env.limits.max_satellite_rotation = deg2rad(10);  % Max allowed θ_sat [rad]
env.limits.max_acceleration = 50;                 % Max joint acceleration [rad/s²]
env.limits.max_velocity = 10;                     % Max joint velocity [rad/s]

%% Display Summary
fprintf('✓ Environment 2D configured:\n');
fprintf('  Time step: %.3f s (%.0f Hz)\n', env.sim.dt, 1/env.sim.dt);
fprintf('  Duration: %.1f s\n', env.sim.T_final);
fprintf('  Contact model: %s (%d points)\n', ...
        env.contact.type, env.contact.n_contacts);
fprintf('  Momentum conservation: %s\n', ...
        char(string(env.momentum.enforce_conservation)));

end
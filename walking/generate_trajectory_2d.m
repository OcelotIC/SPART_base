function traj = generate_trajectory_2d(env, type)
% GENERATE_TRAJECTORY_2D - Creates desired trajectory for walker torso
%
% Inputs:
%   env  - Environment parameters
%   type - Trajectory type: 'line', 'circle', 'static'
%
% Outputs:
%   traj - Structure with fields:
%          .t        - Time vector [N×1]
%          .pos      - Position [3×N]: [x; y; z]
%          .vel      - Velocity [3×N]
%          .acc      - Acceleration [3×N]
%          .q_joints - Desired joint angles [6×N] (from IK)

if nargin < 2
    type = 'line';
end

%% Time vector
dt = env.sim.dt;
T_final = env.sim.T_final;
t = 0:dt:T_final;
N = length(t);

%% Generate trajectory based on type
switch lower(type)
    case 'line'
        % Straight line: move +2m in X direction over duration
        x_start = 0;
        x_end = 2.0;
        z_height = 0.8;  % Constant height
        
        % Use smooth polynomial (5th order) for zero initial/final vel/acc
        s = smooth_step(t / T_final);  % Normalized time [0,1]
        
        x = x_start + (x_end - x_start) * s;
        z = z_height * ones(size(t));
        
        % Derivatives
        ds_dt = smooth_step_derivative(t / T_final) / T_final;
        d2s_dt2 = smooth_step_second_derivative(t / T_final) / T_final^2;
        
        vx = (x_end - x_start) * ds_dt;
        vz = zeros(size(t));
        
        ax = (x_end - x_start) * d2s_dt2;
        az = zeros(size(t));
        
    case 'circle'
        % Circular arc
        radius = 1.0;
        z_center = 0.8;
        omega = 2*pi / T_final;  % One full circle
        
        x = radius * sin(omega * t);
        z = z_center + radius * (1 - cos(omega * t));
        
        vx = radius * omega * cos(omega * t);
        vz = radius * omega * sin(omega * t);
        
        ax = -radius * omega^2 * sin(omega * t);
        az = radius * omega^2 * cos(omega * t);
        
    case 'static'
        % Stationary (for testing dynamics)
        x = zeros(size(t));
        z = 0.8 * ones(size(t));
        vx = zeros(size(t));
        vz = zeros(size(t));
        ax = zeros(size(t));
        az = zeros(size(t));
        
    otherwise
        error('Unknown trajectory type: %s', type);
end

%% Assemble trajectory
traj.t = t;
traj.pos = [x; zeros(size(t)); z];  % Y = 0 (2D)
traj.vel = [vx; zeros(size(t)); vz];
traj.acc = [ax; zeros(size(t)); az];

%% Inverse Kinematics (simplified)
% For each time step, compute desired joint angles
% Assuming fixed contact points at ±1.5m on satellite
traj.q_joints = zeros(6, N);
for i = 1:N
    % Simplified IK: keep arms in "W" configuration
    % This is just for initialization; real IK would solve for contact constraints
    traj.q_joints(:, i) = [
        deg2rad(30);   % arm1_joint1
        deg2rad(-60);  % arm1_joint2
        deg2rad(30);   % arm1_joint3
        deg2rad(-30);  % arm2_joint1
        deg2rad(-60);  % arm2_joint2
        deg2rad(-30)   % arm2_joint3
    ];
end

fprintf('✓ Trajectory generated: %s, %.1f s, %d points\n', type, T_final, N);

end

%% Helper Functions for Smooth Trajectories
function s = smooth_step(t)
    % 5th order polynomial: s(0)=0, s(1)=1, s'(0)=s'(1)=0, s''(0)=s''(1)=0
    s = 6*t.^5 - 15*t.^4 + 10*t.^3;
end

function ds = smooth_step_derivative(t)
    ds = 30*t.^4 - 60*t.^3 + 30*t.^2;
end

function d2s = smooth_step_second_derivative(t)
    d2s = 120*t.^3 - 180*t.^2 + 60*t;
end
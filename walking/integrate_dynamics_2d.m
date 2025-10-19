function state_next = integrate_dynamics_2d(state, tau_joints, sys, env)
% INTEGRATE_DYNAMICS_2D - Time integration of constrained dynamics
%
% Solves the augmented system:
%   [M,  J^T] [q̈  ] = [τ - C·q̇]
%   [J,  0  ] [λ  ]   [-γ_stab ]
%
% where γ_stab = Baumgarte stabilization term
%
% Inputs:
%   state      - Current state
%   tau_joints - Joint torques [6×1]
%   sys        - System dynamics (from compute_coupled_dynamics_2d)
%   env        - Environment parameters
%
% Outputs:
%   state_next - Integrated state at t + dt

dt = env.sim.dt;

%% 1. EXTRACT CURRENT GENERALIZED VELOCITIES
% System: q = [x_w, z_w, θ_w, qm_w(6), x_s, z_s, θ_s]
%         q̇ = [ẋ_w, ż_w, θ̇_w, q̇m_w(6), ẋ_s, ż_s, θ̇_s]

% Walker velocities (2D projection)
u0_walker_2d = [state.u0_walker(2); state.u0_walker(4); state.u0_walker(6)];  % [ωy; vx; vz]
um_walker = state.um_walker;  % [6×1]

% Satellite velocities (2D projection)
u0_sat_2d = [state.u0_sat(2); state.u0_sat(4); state.u0_sat(6)];  % [ωy; vx; vz]

% Full generalized velocity [12×1]
q_dot = [u0_walker_2d; um_walker; u0_sat_2d];

%% 2. GENERALIZED FORCES
% Walker: tau_base = 0 (floating), tau_joints from controller
% Satellite: no actuation (free-floating), RW provides internal torque
tau_full = [zeros(3,1); tau_joints; zeros(3,1)];  % [12×1]

%% 3. BAUMGARTE STABILIZATION
% To prevent constraint drift: γ = -2α·J·q̇ - β²·Φ
alpha_stab = env.contact.alpha_stab;
beta_stab = env.contact.beta_stab;

J = sys.J_constraint;  % [6×12]
Phi = sys.Phi;         % [6×1]

gamma_stab = -2 * alpha_stab * (J * q_dot) - beta_stab^2 * Phi;

%% 4. AUGMENTED SYSTEM
% [M,  J^T] [q̈] = [τ - C·q̇    ]
% [J,  0  ] [λ]   [-γ_stab    ]

M = sys.M_full;  % [12×12]
C = sys.C_full;  % [12×12]
n_q = 12;
n_c = 6;  % Number of constraints

% Right-hand side
rhs_dynamics = tau_full - C * q_dot;
rhs_constraint = -gamma_stab;

% Augmented matrices
M_aug = [M, J'; J, zeros(n_c, n_c)];  % [18×18]
rhs_aug = [rhs_dynamics; rhs_constraint];  % [18×1]

%% 5. SOLVE LINEAR SYSTEM
% Regularization for numerical stability
cond_M_aug = cond(M_aug);
if cond_M_aug > 1e10
    rho = 1e-8;  % Proximal regularization
    M_aug = M_aug + rho * eye(size(M_aug));
    fprintf('  Warning: Ill-conditioned augmented matrix (cond=%.2e), regularizing\n', cond_M_aug);
end

try
    sol = M_aug \ rhs_aug;
catch
    warning('Direct solve failed, using pseudoinverse');
    sol = pinv(M_aug, 1e-8) * rhs_aug;
end

% Extract accelerations and constraint forces
q_ddot = sol(1:n_q);      % [12×1]
lambda = sol(n_q+1:end);  % [6×1]

% Safety clipping
max_accel = env.limits.max_acceleration;
if norm(q_ddot) > max_accel
    q_ddot = q_ddot / norm(q_ddot) * max_accel;
    fprintf('  Warning: Acceleration clipped to %.1f\n', max_accel);
end

%% 6. TIME INTEGRATION (Semi-implicit Euler)
% Update velocities
q_dot_next = q_dot + q_ddot * dt;

% Clip velocities
max_vel = env.limits.max_velocity;
for i = 1:length(q_dot_next)
    if abs(q_dot_next(i)) > max_vel
        q_dot_next(i) = sign(q_dot_next(i)) * max_vel;
    end
end

% Update positions
% Walker base
state_next = state;
state_next.u0_walker(2) = q_dot_next(1);  % ωy
state_next.u0_walker(4) = q_dot_next(2);  % vx
state_next.u0_walker(6) = q_dot_next(3);  % vz

state_next.r0_walker(1) = state.r0_walker(1) + q_dot_next(2) * dt;  % x
state_next.r0_walker(3) = state.r0_walker(3) + q_dot_next(3) * dt;  % z

% Walker orientation (integrate quaternion)
omega_walker = [0; q_dot_next(1); 0];  % Only Y-rotation
if norm(omega_walker) > 1e-10
    % Use SPART quaternion integration if available
    % Otherwise simple Euler: q_next = q + 0.5*dt*Ω*q
    omega_quat = [omega_walker; 0];  % [ωx; ωy; ωz; 0]
    Omega = [0, -omega_walker(1), -omega_walker(2), -omega_walker(3);
             omega_walker(1), 0, omega_walker(3), -omega_walker(2);
             omega_walker(2), -omega_walker(3), 0, omega_walker(1);
             omega_walker(3), omega_walker(2), -omega_walker(1), 0];
    dq = 0.5 * dt * Omega * state.quat_walker;
    state_next.quat_walker = state.quat_walker + dq;
    state_next.quat_walker = state_next.quat_walker / norm(state_next.quat_walker);  % Normalize
else
    state_next.quat_walker = state.quat_walker;
end

% Walker joints
state_next.um_walker = q_dot_next(4:9);
state_next.qm_walker = state.qm_walker + state_next.um_walker * dt;

% Satellite base
state_next.u0_sat(2) = q_dot_next(10);  % ωy
state_next.u0_sat(4) = q_dot_next(11);  % vx
state_next.u0_sat(6) = q_dot_next(12);  % vz

state_next.r0_sat(1) = state.r0_sat(1) + q_dot_next(11) * dt;  % x
state_next.r0_sat(3) = state.r0_sat(3) + q_dot_next(12) * dt;  % z

% Satellite orientation
omega_sat = [0; q_dot_next(10); 0];
if norm(omega_sat) > 1e-10
    Omega_sat = [0, -omega_sat(1), -omega_sat(2), -omega_sat(3);
                 omega_sat(1), 0, omega_sat(3), -omega_sat(2);
                 omega_sat(2), -omega_sat(3), 0, omega_sat(1);
                 omega_sat(3), omega_sat(2), -omega_sat(1), 0];
    dq_sat = 0.5 * dt * Omega_sat * state.quat_sat;
    state_next.quat_sat = state.quat_sat + dq_sat;
    state_next.quat_sat = state_next.quat_sat / norm(state_next.quat_sat);
else
    state_next.quat_sat = state.quat_sat;
end

%% 7. REACTION WHEEL DYNAMICS
% Conservation of angular momentum: ḣ_RW = -ḣ_robot - ḣ_sat
% In 2D, only Z-component matters (perpendicular to plane)

% Compute momentum rate for walker+sat
delta_h_system = sys.h_total(3) - state.h_total_initial(3);  % Should be ≈ 0

% RW compensates (with saturation)
h_RW_desired = -delta_h_system;

% Apply RW torque limit
tau_RW_max = env.satellite_params.RW.tau_max;
tau_RW = max(-tau_RW_max, min(tau_RW_max, -q_ddot(10) * sys.M_full(10,10)));  % Simplified

% Integrate RW momentum
state_next.h_RW = state.h_RW + tau_RW * dt;

% Saturate RW momentum
h_RW_max = env.satellite_params.RW.h_max;
if abs(state_next.h_RW) > h_RW_max
    state_next.h_RW = sign(state_next.h_RW) * h_RW_max;
    fprintf('  Warning: RW saturated at %.2f N·m·s\n', h_RW_max);
end

%% 8. STORE DIAGNOSTICS
state_next.lambda = lambda;  % Constraint forces
state_next.constraint_violation = norm(sys.Phi);
state_next.h_total_initial = state.h_total_initial;  % Preserve initial momentum

end
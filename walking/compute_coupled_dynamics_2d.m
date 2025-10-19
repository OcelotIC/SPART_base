function sys = compute_coupled_dynamics_2d(state, robot_walker, robot_sat, satellite_params, env)
% COMPUTE_COUPLED_DYNAMICS_2D - Computes coupled spacecraft-robot dynamics
%
% Assembles the augmented mass matrix, Coriolis matrix, and constraint Jacobian
% for the complete system: walker (9 DoF) + satellite (3 DoF) = 12 DoF total
%
% Inputs:
%   state            - Current state structure
%   robot_walker     - Walker SPART model
%   robot_sat        - Satellite SPART model
%   satellite_params - Satellite parameters
%   env              - Environment parameters
%
% Outputs:
%   sys - Structure containing:
%         .M_full       - Mass matrix [12×12]
%         .C_full       - Coriolis matrix [12×12]
%         .J_constraint - Constraint Jacobian [6×12]
%         .Phi          - Constraint violation [6×1]
%         .h_total      - Total angular momentum [3×1]
%         .kinematic    - Forward kinematics data

%% 1. WALKER DYNAMICS (SPART)
% Extract 2D state from full 3D representation
% SPART quaternion: [qx; qy; qz; qw] where qw is scalar
R0_walker = quat_DCM(state.quat_walker);  % Use SPART native function
r0_walker = state.r0_walker;
qm_walker = state.qm_walker;
u0_walker = state.u0_walker;  % [ωx; ωy; ωz; vx; vy; vz]
um_walker = state.um_walker;

% Kinematics
[RJ_w, RL_w, rJ_w, rL_w, e_w, g_w] = Kinematics(R0_walker, r0_walker, qm_walker, robot_walker);

% Differential kinematics
[Bij_w, Bi0_w, P0_w, pm_w] = DiffKinematics(R0_walker, r0_walker, rL_w, e_w, g_w, robot_walker);

% Velocities
[t0_w, tm_w] = Velocities(Bij_w, Bi0_w, P0_w, pm_w, u0_walker, um_walker, robot_walker);

% Inertias in inertial frame
[I0_w, Im_w] = I_I(R0_walker, RL_w, robot_walker);

% Mass composite body
[M0_tilde_w, Mm_tilde_w] = MCB(I0_w, Im_w, Bij_w, Bi0_w, robot_walker);

% Generalized inertia matrix
[H0_w, H0m_w, Hm_w] = GIM(M0_tilde_w, Mm_tilde_w, Bij_w, Bi0_w, P0_w, pm_w, robot_walker);

% Coriolis matrix
[C0_w, C0m_w, Cm0_w, Cm_w] = CIM(t0_w, tm_w, I0_w, Im_w, M0_tilde_w, Mm_tilde_w, ...
                                  Bij_w, Bi0_w, P0_w, pm_w, robot_walker);

% Walker full matrices [9×9]
M_walker = [H0_w, H0m_w; H0m_w', Hm_w];  % [6+6 × 6+6]
C_walker = [C0_w, C0m_w; Cm0_w, Cm_w];

%% 2. SATELLITE DYNAMICS (Rigid Body)
R0_sat = quat_DCM(state.quat_sat);
r0_sat = state.r0_sat;
u0_sat = state.u0_sat;  % [ωx; ωy; ωz; vx; vy; vz]

% Mass matrix (rigid body, 6×6)
m_sat = robot_sat.base_link.mass;
I_sat = robot_sat.base_link.inertia(1:3, 1:3);  % Rotational inertia
M_sat = [I_sat, zeros(3,3); zeros(3,3), m_sat * eye(3)];

% Coriolis matrix (rigid body)
omega_sat = u0_sat(1:3);
v_sat = u0_sat(4:6);
omega_skew = skew(omega_sat);
C_sat = [omega_skew * I_sat, zeros(3,3); 
         zeros(3,3), m_sat * omega_skew];

%% 3. EXTRACT 2D COMPONENTS
% System has 12 DoF in 2D: [x_w, z_w, θ_w, qm_w(6), x_s, z_s, θ_s]
% But SPART works in 6D. We need to project to 2D.

% Define projection matrix: 6D → 3D (2D + rotation around Y)
% From u0 = [ωx; ωy; ωz; vx; vy; vz] to [θ̇; ẋ; ż]
% Map: θ ↔ ωy(2), x ↔ vx(4), z ↔ vz(6)
P_2d = zeros(3, 6);
P_2d(1, 2) = 1;  % θ̇ = ωy
P_2d(2, 4) = 1;  % ẋ = vx
P_2d(3, 6) = 1;  % ż = vz

% Walker 2D matrices [9×9]
M_walker_2d = zeros(9, 9);
M_walker_2d(1:3, 1:3) = P_2d * H0_w * P_2d';  % Base-base
M_walker_2d(1:3, 4:9) = P_2d * H0m_w;         % Base-joints
M_walker_2d(4:9, 1:3) = H0m_w' * P_2d';       % Joints-base
M_walker_2d(4:9, 4:9) = Hm_w;                 % Joints-joints

C_walker_2d = zeros(9, 9);
C_walker_2d(1:3, 1:3) = P_2d * C0_w * P_2d';
C_walker_2d(1:3, 4:9) = P_2d * C0m_w;
C_walker_2d(4:9, 1:3) = Cm0_w * P_2d';
C_walker_2d(4:9, 4:9) = Cm_w;

% Satellite 2D matrices [3×3]
M_sat_2d = P_2d * M_sat * P_2d';
C_sat_2d = P_2d * C_sat * P_2d';

%% 4. COUPLED SYSTEM MATRICES [12×12]
sys.M_full = zeros(12, 12);
sys.M_full(1:9, 1:9) = M_walker_2d;      % Walker block
sys.M_full(10:12, 10:12) = M_sat_2d;     % Satellite block
% No coupling in mass matrix (coupled through constraints)

sys.C_full = zeros(12, 12);
sys.C_full(1:9, 1:9) = C_walker_2d;
sys.C_full(10:12, 10:12) = C_sat_2d;

%% 5. CONTACT CONSTRAINTS
[Phi, J_constraint, r_ee] = compute_contact_constraints_2d(...
    R0_walker, r0_walker, qm_walker, robot_walker, ...
    R0_sat, r0_sat, satellite_params);

sys.Phi = Phi;
sys.J_constraint = J_constraint;

%% 6. MOMENTUM (for conservation check)
% Total angular momentum (Z-component only in 2D)
% h_total = I_walker * ω_walker + I_sat * ω_sat + h_RW

% Walker angular momentum (around Y-axis for 2D)
h_walker_y = (P_2d * H0_w * P_2d') * [state.u0_walker(2); state.u0_walker(4); state.u0_walker(6)] + ...
             (P_2d * H0m_w) * state.um_walker;
h_walker = [0; h_walker_y(1); 0];  % Only Y-component matters in 2D

% Satellite angular momentum
h_sat_y = M_sat_2d(1,1) * state.u0_sat(2);  % I_sat * ω_sat_y
h_sat = [0; h_sat_y; 0];

% Reaction wheel momentum
h_RW = [0; 0; state.h_RW];  % Z-component (perpendicular to 2D plane)

% Total
sys.h_total = h_walker + h_sat + h_RW;

%% 7. KINEMATICS DATA (for visualization)
sys.kinematic.walker.RL = RL_w;
sys.kinematic.walker.rL = rL_w;
sys.kinematic.walker.r_ee = r_ee;
sys.kinematic.satellite.r0 = r0_sat;
sys.kinematic.satellite.R0 = R0_sat;

end

%% Helper Functions
function S = skew(v)
    % Skew-symmetric matrix from 3D vector
    S = [0, -v(3), v(2);
         v(3), 0, -v(1);
         -v(2), v(1), 0];
end
function [Phi, J_constraint, r_ee] = compute_contact_constraints_2d(...
    R0_walker, r0_walker, qm_walker, robot_walker, ...
    R0_sat, r0_sat, satellite_params)
% COMPUTE_CONTACT_CONSTRAINTS_2D - Bilateral contact constraints for docking
%
% Computes position and orientation constraints for end-effectors
% docked to satellite surface (3 constraints per contact × 2 contacts = 6 total)
%
% Inputs:
%   R0_walker, r0_walker, qm_walker - Walker state
%   robot_walker                    - Walker SPART model
%   R0_sat, r0_sat                  - Satellite state
%   satellite_params                - Satellite parameters (contact points)
%
% Outputs:
%   Phi          - Constraint violation [6×1]: [Φ1_x; Φ1_z; Φ1_θ; Φ2_x; Φ2_z; Φ2_θ]
%   J_constraint - Constraint Jacobian [6×15]: ∂Φ/∂q_system
%   r_ee         - End-effector positions [3×2] for visualization

%% 1. Kinematics - Walker
% Compute forward kinematics for all links
[RJ, RL, rJ, rL, e, g] = Kinematics(R0_walker, r0_walker, qm_walker, robot_walker);

% End-effector positions (links 3 and 6)
% Link 3: end of arm 1 (right)
% Link 6: end of arm 2 (left)
r_ee1_inertial = rL(:, 3);  % End-effector 1 position [3×1]
R_ee1_inertial = RL(:, :, 3);  % End-effector 1 orientation

r_ee2_inertial = rL(:, 6);  % End-effector 2 position [3×1]
R_ee2_inertial = RL(:, :, 6);

% Extract orientation angle (2D: rotation around Y-axis)
% For 2D in XZ plane, θ = atan2(R(3,1), R(1,1))
theta_ee1 = atan2(R_ee1_inertial(3,1), R_ee1_inertial(1,1));
theta_ee2 = atan2(R_ee2_inertial(3,1), R_ee2_inertial(1,1));

%% 2. Kinematics - Satellite Contact Points
% Contact point 1 (right, +X)
r_contact1_local = satellite_params.contact_points(1).position_local;
r_contact1_inertial = r0_sat + R0_sat * r_contact1_local;
theta_sat = atan2(R0_sat(3,1), R0_sat(1,1));
theta_contact1 = theta_sat + satellite_params.contact_points(1).orientation_local;

% Contact point 2 (left, -X)
r_contact2_local = satellite_params.contact_points(2).position_local;
r_contact2_inertial = r0_sat + R0_sat * r_contact2_local;
theta_contact2 = theta_sat + satellite_params.contact_points(2).orientation_local;

%% 3. Constraint Violation Φ(q)
% Each contact: 3 constraints (x, z, θ)
% Φ = 0 means perfect docking

Phi = zeros(6, 1);

% Contact 1 constraints
Phi(1) = r_ee1_inertial(1) - r_contact1_inertial(1);  % x position
Phi(2) = r_ee1_inertial(3) - r_contact1_inertial(3);  % z position
Phi(3) = wrap_angle(theta_ee1 - theta_contact1);      % θ orientation

% Contact 2 constraints
Phi(4) = r_ee2_inertial(1) - r_contact2_inertial(1);
Phi(5) = r_ee2_inertial(3) - r_contact2_inertial(3);
Phi(6) = wrap_angle(theta_ee2 - theta_contact2);

%% 4. Constraint Jacobian J = ∂Φ/∂q
% System coordinates: q_system = [x_walker; z_walker; θ_walker; qm_walker(6); 
%                                  x_sat; z_sat; θ_sat]
% Total: 3 + 6 + 3 = 12 DoF (in 2D, y components ignored)

% Use SPART's Jacobian for walker end-effectors
[J0_ee1, Jm_ee1] = Jacob(r_ee1_inertial, r0_walker, rL, ...
                         zeros(6,6), zeros(6,robot_walker.n_q), ...
                         3, robot_walker);  % Link 3

[J0_ee2, Jm_ee2] = Jacob(r_ee2_inertial, r0_walker, rL, ...
                         zeros(6,6), zeros(6,robot_walker.n_q), ...
                         6, robot_walker);  % Link 6

% Extract 2D components (x, z from position; ωy for rotation)
% J0 is [6×6]: [ω(3); v(3)]
% Jm is [6×n_q]

% For position (x, z):
J0_ee1_pos_xz = J0_ee1([4,6], :);  % vx, vz from translational part
Jm_ee1_pos_xz = Jm_ee1([4,6], :);

J0_ee2_pos_xz = J0_ee2([4,6], :);
Jm_ee2_pos_xz = Jm_ee2([4,6], :);

% For orientation (θ around Y):
J0_ee1_orient = J0_ee1(2, :);  % ωy
Jm_ee1_orient = Jm_ee1(2, :);

J0_ee2_orient = J0_ee2(2, :);
Jm_ee2_orient = Jm_ee2(2, :);

% Satellite Jacobian (simple rigid body)
% ∂r_contact/∂q_sat = [I, 0; 0, I; lever_arm_z, lever_arm_x, 1]
theta_sat_current = atan2(R0_sat(3,1), R0_sat(1,1));

J_sat_contact1 = zeros(3, 3);  % [∂x/∂q_sat, ∂z/∂q_sat, ∂θ/∂q_sat]
J_sat_contact1(1, 1) = 1;  % x_contact depends on x_sat
J_sat_contact1(2, 2) = 1;  % z_contact depends on z_sat
J_sat_contact1(1, 3) = -r_contact1_local(3) * sin(theta_sat_current) - ...
                        r_contact1_local(1) * cos(theta_sat_current);
J_sat_contact1(2, 3) = r_contact1_local(3) * cos(theta_sat_current) - ...
                       r_contact1_local(1) * sin(theta_sat_current);
J_sat_contact1(3, 3) = 1;  % θ_contact depends on θ_sat

J_sat_contact2 = zeros(3, 3);
J_sat_contact2(1, 1) = 1;
J_sat_contact2(2, 2) = 1;
J_sat_contact2(1, 3) = -r_contact2_local(3) * sin(theta_sat_current) - ...
                        r_contact2_local(1) * cos(theta_sat_current);
J_sat_contact2(2, 3) = r_contact2_local(3) * cos(theta_sat_current) - ...
                       r_contact2_local(1) * sin(theta_sat_current);
J_sat_contact2(3, 3) = 1;

% Assemble full Jacobian [6×12]
% Order: [walker_base(3), walker_joints(6), sat_base(3)]
% Need to extract 2D components from J0 (which is 6D)

% Simplified: extract only relevant 2D DoF mapping
% Walker base: [x, z, θ] → extract columns for vx, vz, ωy
% From u0 = [ωx; ωy; ωz; vx; vy; vz], we need [ωy(2); vx(4); vz(6)]
J0_2d_map = [2; 4; 6];  % Indices for θ, x, z

J_constraint = zeros(6, 12);

% Contact 1: rows 1-3
J_constraint(1:2, 1:3) = J0_ee1_pos_xz(:, J0_2d_map);  % Walker base contribution
J_constraint(1:2, 4:9) = Jm_ee1_pos_xz;                % Walker joints
J_constraint(1:2, 10:12) = -J_sat_contact1(1:2, :);    % Satellite (negative)

J_constraint(3, 1:3) = J0_ee1_orient(J0_2d_map);
J_constraint(3, 4:9) = Jm_ee1_orient;
J_constraint(3, 10:12) = -J_sat_contact1(3, :);

% Contact 2: rows 4-6
J_constraint(4:5, 1:3) = J0_ee2_pos_xz(:, J0_2d_map);
J_constraint(4:5, 4:9) = Jm_ee2_pos_xz;
J_constraint(4:5, 10:12) = -J_sat_contact2(1:2, :);

J_constraint(6, 1:3) = J0_ee2_orient(J0_2d_map);
J_constraint(6, 4:9) = Jm_ee2_orient;
J_constraint(6, 10:12) = -J_sat_contact2(3, :);

%% 5. Output end-effector positions for visualization
r_ee = [r_ee1_inertial, r_ee2_inertial];

end

%% Helper Function
function theta = wrap_angle(theta)
    % Wrap angle to [-π, π]
    theta = mod(theta + pi, 2*pi) - pi;
end
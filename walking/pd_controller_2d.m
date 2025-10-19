function tau_joints = pd_controller_2d(state, q_desired, qd_desired, env)
% PD_CONTROLLER_2D - Simple PD controller for joint space control
%
% Computes joint torques: τ = Kp·(q_des - q) + Kd·(q̇_des - q̇)
%
% Inputs:
%   state      - Current state
%   q_desired  - Desired joint positions [6×1]
%   qd_desired - Desired joint velocities [6×1]
%   env        - Environment parameters
%
% Outputs:
%   tau_joints - Joint torques [6×1]

%% Extract current joint state
q_current = state.qm_walker;
qd_current = state.um_walker;

%% PD gains from environment
Kp = env.controller.Kp;
Kd = env.controller.Kd;

%% Compute torque
tau_joints = Kp * (q_desired - q_current) + Kd * (qd_desired - qd_current);

%% Torque limits (optional)
tau_max = 100;  % N·m per joint
for i = 1:length(tau_joints)
    if abs(tau_joints(i)) > tau_max
        tau_joints(i) = sign(tau_joints(i)) * tau_max;
    end
end

end
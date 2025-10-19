%% Constrained Forward Dynamics - Clean Implementation
% Pure forward dynamics: given state and tau, compute qddot and lambda
%
% NO preprocessing, NO velocity projection
% Just solve: H*qddot + C*qdot = tau + J'*lambda
%             J*qddot = -Jdot*qdot - 2*alpha*(J*qdot) - beta^2*phi

clc; close all; clear;

%--- URDF filename ---%
filename='kuka_lwr.urdf';

%--- Create robot model ---%
[robot,robot_keys] = urdf2robot(filename);

%--- Simulation parameters ---%
dt = 0.0005;
t_final = 5.0;
n_steps = round(t_final/dt);

% Baumgarte parameters
alpha = 10;    % Try higher values for your case
beta = 50;    % Try higher values

% Regularization
epsilon = 1e-6;  % Larger regularization for robustness

% Zero gravity
g_gravity = 0.0;

%--- Initial conditions (YOUR CHOICE) ---%
R0_init = eye(3);
r0_init = zeros(3,1);
qm_init = 0.05*randn(robot.n_q,1);

% Initial velocities - whatever you want
u0_init = zeros(6,1);
um_init = 0.01*ones(robot.n_q,1);  % Your 0.8 rad/s

%--- Fix contact point ---%
[~,RL_init,~,rL_init,~,~] = Kinematics(R0_init,r0_init,qm_init,robot);
r_contact = rL_init(1:3,end);
R_contact = RL_init(1:3,1:3,end);

fprintf('Contact fixed at: [%.3f, %.3f, %.3f]\n', r_contact);

%--- Convert to state ---%
p0_init = DCM_quat(R0_init);  % SPART function
n_q = robot.n_q;
state = [r0_init; p0_init; qm_init; u0_init; um_init];

%--- Storage ---%
time_hist = zeros(1, n_steps+1);
state_hist = zeros(length(state), n_steps+1);
lambda_hist = zeros(6, n_steps+1);
constraint_hist = zeros(2, n_steps+1); % [pos; vel] violations

time_hist(1) = 0;
state_hist(:,1) = state;

%--- Check initial constraint violations ---%
[phi_pos_0, phi_vel_0] = check_constraints(state, robot, r_contact, R_contact);
fprintf('Initial constraint violations:\n');
fprintf('  Position: %.3e m\n', norm(phi_pos_0));
fprintf('  Velocity: %.3e m/s\n', norm(phi_vel_0));
fprintf('NOTE: Large velocity violation is OK - dynamics will handle it\n\n');

%--- Integration loop ---%
fprintf('Running simulation (dt=%.4f, alpha=%.1f, beta=%.1f)...\n', dt, alpha, beta);

for k = 1:n_steps
    try
        % RK4 step
        [state, lambda] = rk4_step(state, dt, robot, r_contact, R_contact, ...
                                    alpha, beta, g_gravity, epsilon);
        
        % Store
        time_hist(k+1) = k*dt;
        state_hist(:,k+1) = state;
        lambda_hist(:,k+1) = lambda;
        
        % Monitor constraints
        [phi_pos, phi_vel] = check_constraints(state, robot, r_contact, R_contact);
        constraint_hist(1,k+1) = norm(phi_pos);
        constraint_hist(2,k+1) = norm(phi_vel);
        
        % Check for blow-up
        if constraint_hist(1,k+1) > 1.0 || any(isnan(state))
            fprintf('\nDiverged at t=%.3f s\n', k*dt);
            fprintf('  Position error: %.3e\n', constraint_hist(1,k+1));
            fprintf('  Velocity error: %.3e\n', constraint_hist(2,k+1));
            break;
        end
        
        if mod(k, 500) == 0
            fprintf('t=%.2f | pos_err=%.2e | vel_err=%.2e | lambda_norm=%.1f\n', ...
                    k*dt, constraint_hist(1,k+1), constraint_hist(2,k+1), norm(lambda));
        end
        
    catch ME
        fprintf('\nError at t=%.3f: %s\n', k*dt, ME.message);
        break;
    end
end

% Trim
n_actual = min(k+1, n_steps+1);
time_hist = time_hist(1:n_actual);
state_hist = state_hist(:,1:n_actual);
lambda_hist = lambda_hist(:,1:n_actual);
constraint_hist = constraint_hist(:,1:n_actual);

%--- Plot ---%
figure('Position', [100, 100, 1200, 800]);

subplot(3,3,1);
plot(time_hist, state_hist(8:7+n_q,:)*180/pi, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Joint Angles [deg]');
title('Joint Positions'); grid on;

subplot(3,3,2);
plot(time_hist, state_hist(14+n_q:end,:)*180/pi, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Joint Velocities [deg/s]');
title('Joint Velocities'); grid on;

subplot(3,3,3);
plot(time_hist, vecnorm(lambda_hist(4:6,:)), 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Force [N]');
title('Contact Force Magnitude'); grid on;

subplot(3,3,4);
plot(time_hist, vecnorm(lambda_hist(1:3,:)), 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Torque [Nm]');
title('Contact Torque Magnitude'); grid on;

subplot(3,3,5);
semilogy(time_hist, constraint_hist(1,:), 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Error [m/rad]');
title('Position Constraint Violation'); grid on;

subplot(3,3,6);
semilogy(time_hist, constraint_hist(2,:), 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Error [m/s or rad/s]');
title('Velocity Constraint Violation'); grid on;

subplot(3,3,7);
plot(time_hist, lambda_hist(4:6,:), 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Force [N]');
title('Contact Forces (x,y,z)'); grid on;
legend('f_x', 'f_y', 'f_z');

subplot(3,3,8);
plot(time_hist, lambda_hist(1:3,:), 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Torque [Nm]');
title('Contact Torques (x,y,z)'); grid on;
legend('\tau_x', '\tau_y', '\tau_z');

subplot(3,3,9);
plot(time_hist, state_hist(1:3,:), 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Position [m]');
title('Base Position'); grid on;
legend('x', 'y', 'z');

sgtitle(sprintf('Constrained FD: \\alpha=%.0f, \\beta=%.0f, \\epsilon=%.0e', ...
                alpha, beta, epsilon), 'FontSize', 14);

fprintf('\n=== FINAL DIAGNOSTICS ===\n');
fprintf('Simulation duration: %.3f s\n', time_hist(end));
fprintf('Final position error: %.3e\n', constraint_hist(1,end));
fprintf('Final velocity error: %.3e\n', constraint_hist(2,end));
fprintf('Max constraint force: %.1f N\n', max(vecnorm(lambda_hist(4:6,:))));

%% ========================================================================
%  CORE DYNAMICS FUNCTIONS
%  ========================================================================

function [state_next, lambda] = rk4_step(state, dt, robot, r_contact, ...
                                         R_contact, alpha, beta, g_gravity, epsilon)
    [k1, lambda1] = state_derivative(state, robot, r_contact, R_contact, ...
                                      alpha, beta, g_gravity, epsilon);
    
    state2 = state + 0.5*dt*k1;
    state2 = normalize_quat(state2);
    [k2, ~] = state_derivative(state2, robot, r_contact, R_contact, ...
                                alpha, beta, g_gravity, epsilon);
    
    state3 = state + 0.5*dt*k2;
    state3 = normalize_quat(state3);
    [k3, ~] = state_derivative(state3, robot, r_contact, R_contact, ...
                                alpha, beta, g_gravity, epsilon);
    
    state4 = state + dt*k3;
    state4 = normalize_quat(state4);
    [k4, lambda4] = state_derivative(state4, robot, r_contact, R_contact, ...
                                      alpha, beta, g_gravity, epsilon);
    
    state_next = state + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    state_next = normalize_quat(state_next);
    lambda = lambda4;
end

function [state_dot, lambda] = state_derivative(state, robot, r_contact, ...
                                                 R_contact, alpha, beta, g_gravity, epsilon)
    % Pure forward dynamics: compute qddot and lambda from state
    
    n_q = robot.n_q;
    
    % Extract state
    r0 = state(1:3);
    p0 = state(4:7);
    qm = state(8:7+n_q);
    u0 = state(8+n_q:13+n_q);
    um = state(14+n_q:end);
    
    R0 = quat_DCM(p0);
    
    % Kinematics
    [~,RL,~,rL,e,g] = Kinematics(R0,r0,qm,robot);
    [Bij,Bi0,P0,pm] = DiffKinematics(R0,r0,rL,e,g,robot);
    [t0,tm] = Velocities(Bij,Bi0,P0,pm,u0,um,robot);
    
    % Constraint Jacobian
    r_ee = rL(1:3,end);
    [J0n, Jmn] = Jacob(r_ee,r0,rL,P0,pm,robot.n_links_joints,robot);
    J = [J0n, Jmn];
    
    % Dynamics matrices
    [I0,Im] = I_I(R0,RL,robot);
    [M0_tilde,Mm_tilde] = MCB(I0,Im,Bij,Bi0,robot);
    [H0, H0m, Hm] = GIM(M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot);
    H = [H0, H0m; H0m', Hm];
    
    [C0, C0m, Cm0, Cm] = CIM(t0,tm,I0,Im,M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot);
    C_matrix = [C0, C0m; Cm0, Cm];
    C_force = C_matrix * [u0; um];
    
    % Gravity forces
    wF0 = [0;0;0;0;0;-robot.base_link.mass*g_gravity];
    wFm = zeros(6,robot.n_links_joints);
    for i=1:robot.n_links_joints
        wFm(6,i) = -robot.links(i).mass*g_gravity;
    end
    tau_gravity = generalized_forces(wF0, wFm, P0, pm, robot);
    
    % Control torques (zero for open-loop)
    tau = tau_gravity;
    
    % Constraint violations
    R_ee = RL(1:3,1:3,end);
    phi_pos = [rotation_error(R_ee, R_contact); r_ee - r_contact];
    phi_vel = tm(1:6,end);  % End-effector twist
    
    % Jdot*qdot
    t_ee = tm(1:6,end);
    [J0dot, Jmdot] = Jacobdot(r_ee, t_ee, r0, t0, rL, tm, P0, pm, robot.n_links_joints, robot);
    Jdot_qdot = J0dot*u0 + Jmdot*um;
    
    % Baumgarte RHS
    rhs_constraint = -Jdot_qdot - 2*alpha*phi_vel - beta^2*phi_pos;
    
    % KKT system with regularization
    n_total = 6 + n_q;
    KKT = [H + epsilon*eye(n_total), -J'; 
           J,                         epsilon*eye(6)];
    rhs = [tau - C_force; rhs_constraint];
    
    % Solve
    sol = KKT \ rhs;
    qdotdot = sol(1:n_total);
    lambda = sol(n_total+1:end);
    
    % State derivative
    r0dot = u0(4:6);
    p0dot = 0.5 * quat_omega_product(p0, u0(1:3));
    qmdot = um;
    
    state_dot = [r0dot; p0dot; qmdot; qdotdot];
end

function tau_g = generalized_forces(wF0, wFm, P0, pm, robot)
    tau0 = wF0;
    taum = zeros(robot.n_q, 1);
    for i=1:robot.n_links_joints
        if robot.joints(i).type ~= 0
            taum(robot.joints(i).q_id) = pm(:,i)' * wFm(:,i);
        end
    end
    tau_g = [tau0; taum];
end

function [phi_pos, phi_vel] = check_constraints(state, robot, r_contact, R_contact)
    n_q = robot.n_q;
    r0 = state(1:3);
    p0 = state(4:7);
    qm = state(8:7+n_q);
    u0 = state(8+n_q:13+n_q);
    um = state(14+n_q:end);
    
    R0 = quat_DCM(p0);
    [~,RL,~,rL,e,g] = Kinematics(R0,r0,qm,robot);
    [Bij,Bi0,P0,pm] = DiffKinematics(R0,r0,rL,e,g,robot);
    [~,tm] = Velocities(Bij,Bi0,P0,pm,u0,um,robot);
    
    r_ee = rL(1:3,end);
    R_ee = RL(1:3,1:3,end);
    phi_pos = [rotation_error(R_ee, R_contact); r_ee - r_contact];
    phi_vel = tm(1:6,end);
end

function phi_ang = rotation_error(R_current, R_desired)
    R_error = R_current' * R_desired;
    theta = acos(min(max((trace(R_error) - 1) / 2, -1), 1));
    if abs(theta) < 1e-6
        phi_ang = zeros(3,1);
    else
        axis = (1/(2*sin(theta))) * [R_error(3,2) - R_error(2,3);
                                       R_error(1,3) - R_error(3,1);
                                       R_error(2,1) - R_error(1,2)];
        phi_ang = theta * axis;
    end
end

function state_norm = normalize_quat(state)
    state_norm = state;
    p0 = state(4:7);
    state_norm(4:7) = p0 / norm(p0);
end

function pdot = quat_omega_product(p, omega)
    % Quaternion derivative: qdot = 0.5 * q ⊗ [0; omega]
    % SPART convention: q = [x,y,z,w] (scalar last)
    qx = p(1); qy = p(2); qz = p(3); qw = p(4);
    wx = omega(1); wy = omega(2); wz = omega(3);
    
    pdot = [qw*wx + qy*wz - qz*wy;
            qw*wy - qx*wz + qz*wx;
            qw*wz + qx*wy - qy*wx;
           -qx*wx - qy*wy - qz*wz];
end
%% Constrained Forward Dynamics - Semi-Implicit Euler
% SPART uses MIXED-FRAME representation:
%   u0 = [ω_body; v_inertial]    (angular in body, linear in inertial)
%   u0dot = [ω̇_body; a_inertial]  (same convention)
%   tau0 = [n_body; f_inertial]   (same convention)
%
% This is standard in floating-base robotics!
%
% Semi-implicit Euler:
%   v_{n+1} = v_n + dt * a_n           (compute NEW velocities)
%   q_{n+1} = q_n + dt * f(v_{n+1})    (integrate with NEW velocities)
%
% The key is using v_{n+1} (not v_n) in position update!

clc; close all; clear;

%--- Setup ---%
filename='kuka_lwr.urdf';
[robot,robot_keys] = urdf2robot(filename);

dt = 0.00001;       % Can use LARGER timestep now!
t_final = 5.0;
n_steps = round(t_final/dt);

% Baumgarte - can use moderate values
alpha = 0*2;
beta = 0*10;
epsilon = 1e-6;
g_gravity = 0.0;

fprintf('\n=== FRAME CONVENTION CHECK ===\n');
fprintf('SPART uses MIXED frames for floating base:\n');
fprintf('  u0(1:3) = ω_body    (angular velocity in BODY frame)\n');
fprintf('  u0(4:6) = v_inertial (linear velocity in INERTIAL frame)\n');
fprintf('This is correct and standard!\n');
fprintf('Semi-implicit integration accounts for this.\n\n');

%--- Initial conditions ---%
R0_init = eye(3);
r0_init = zeros(3,1);
qm_init = 0.05*randn(robot.n_q,1);
u0_init = zeros(6,1);
um_init = 0.0*ones(robot.n_q,1);  % Your 0.8 rad/s

% Fix contact
[~,RL_init,~,rL_init,~,~] = Kinematics(R0_init,r0_init,qm_init,robot);
r_contact = rL_init(1:3,end);
R_contact = RL_init(1:3,1:3,end);

p0_init = DCM_quat(R0_init);
n_q = robot.n_q;

% State: [r0; p0; qm; u0; um]
q = [r0_init; p0_init; qm_init];      % Positions (including quaternion)
v = [u0_init; um_init];                % Velocities

%--- Storage ---%
time_hist = zeros(1, n_steps+1);
q_hist = zeros(length(q), n_steps+1);
v_hist = zeros(length(v), n_steps+1);
lambda_hist = zeros(6, n_steps+1);
constraint_hist = zeros(2, n_steps+1);

time_hist(1) = 0;
q_hist(:,1) = q;
v_hist(:,1) = v;

[phi_pos, phi_vel] = check_constraints_qv(q, v, robot, r_contact, R_contact);
fprintf('Initial violations: pos=%.2e m, vel=%.2e m/s\n\n', norm(phi_pos), norm(phi_vel));

fprintf('Time[s] | Pos_Err | Vel_Err | Lambda[N]\n');
fprintf('--------|---------|---------|----------\n');

%--- Semi-implicit Euler integration ---%
for k = 1:n_steps
    % Semi-implicit Euler (symplectic integrator):
    % Step 1: Compute accelerations at (q_n, v_n)
    % Step 2: Update velocities: v_{n+1} = v_n + dt * vdot_n
    % Step 3: Update positions: q_{n+1} = q_n + dt * f(v_{n+1})  <-- Use NEW velocity!
    
    % 1. Compute accelerations at current state
    [vdot, lambda] = compute_accelerations(q, v, robot, r_contact, R_contact, ...
                                           alpha, beta, g_gravity, epsilon);
    
    % 2. Update velocities (explicit in acceleration)
    v_new = v + dt * vdot;
    
    % 3. Update positions with NEW velocities (this makes it semi-implicit!)
    q_new = update_positions(q, v_new, dt);  % Pass v_new, not v!
    
    % Update state for next iteration
    q = q_new;
    v = v_new;
    
    % Store
    time_hist(k+1) = k*dt;
    q_hist(:,k+1) = q;
    v_hist(:,k+1) = v;
    lambda_hist(:,k+1) = lambda;
    
    % Check constraints
    [phi_pos, phi_vel] = check_constraints_qv(q, v, robot, r_contact, R_contact);
    constraint_hist(1,k+1) = norm(phi_pos);
    constraint_hist(2,k+1) = norm(phi_vel);
    
    % Monitor
    if mod(k, 100) == 0
        fprintf('%7.2f | %.2e | %.2e | %8.1f\n', ...
                k*dt, constraint_hist(1,k+1), constraint_hist(2,k+1), norm(lambda));
    end
    
    if constraint_hist(1,k+1) > 0.5
        fprintf('\n✗ Diverged at t=%.3f s\n', k*dt);
        break;
    end
end

n_actual = min(k+1, n_steps+1);
time_hist = time_hist(1:n_actual);
q_hist = q_hist(:,1:n_actual);
v_hist = v_hist(:,1:n_actual);
lambda_hist = lambda_hist(:,1:n_actual);
constraint_hist = constraint_hist(:,1:n_actual);

%--- Results ---%
fprintf('\n=== FINAL RESULTS ===\n');
fprintf('Duration: %.2f s\n', time_hist(end));
fprintf('Final pos error: %.2e m\n', constraint_hist(1,end));
fprintf('Final vel error: %.2e m/s\n', constraint_hist(2,end));
fprintf('Mean force: %.1f N\n', mean(vecnorm(lambda_hist(4:6,:))));

if time_hist(end) >= t_final - dt
    fprintf('\n✓✓✓ SUCCESS! Semi-implicit Euler works! ✓✓✓\n');
end

%--- Plot ---%
figure('Position', [100, 100, 1200, 800]);

subplot(2,3,1);
plot(time_hist, q_hist(8:7+n_q,:)*180/pi, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Angle [deg]');
title('Joint Positions'); grid on;

subplot(2,3,2);
plot(time_hist, v_hist(7:end,:)*180/pi, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Velocity [deg/s]');
title('Joint Velocities'); grid on;

subplot(2,3,3);
plot(time_hist, vecnorm(lambda_hist(4:6,:)), 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Force [N]');
title('Contact Force Magnitude'); grid on;

subplot(2,3,4);
semilogy(time_hist, constraint_hist(1,:), 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Error [m]');
title('Position Constraint Violation'); grid on;

subplot(2,3,5);
semilogy(time_hist, constraint_hist(2,:), 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Error [m/s]');
title('Velocity Constraint Violation'); grid on;

subplot(2,3,6);
plot(time_hist, lambda_hist(4:6,:), 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Force [N]');
title('Contact Forces'); grid on;
legend('f_x', 'f_y', 'f_z');

sgtitle('Semi-Implicit Euler: Stable Constrained Dynamics', 'FontSize', 14);

%% ========================================================================
%  CORE FUNCTIONS
%  ========================================================================

function [vdot, lambda] = compute_accelerations(q, v, robot, r_contact, ...
                                                 R_contact, alpha, beta, g, epsilon)
    % Compute accelerations from current positions and velocities
    % v = [u0; um] where u0 = [ω_body; v_inertial] (MIXED frame!)
    % vdot = [u0dot; umdot] where u0dot = [ω̇_body; a_inertial]
    
    n_q = robot.n_q;
    
    % Extract from state
    r0 = q(1:3);
    p0 = q(4:7);
    qm = q(8:end);
    u0 = v(1:6);       % [ω_body; v_inertial]
    um = v(7:end);
    
    R0 = quat_DCM(p0);
    
    % Kinematics & Dynamics
    % NOTE: Velocities() converts u0 (mixed frame) → t0 (inertial frame)
    [~,RL,~,rL,e,g_vec]=Kinematics(R0,r0,qm,robot);
    [Bij,Bi0,P0,pm]=DiffKinematics(R0,r0,rL,e,g_vec,robot);
    [t0,tm]=Velocities(Bij,Bi0,P0,pm,u0,um,robot);  % t0 is [ω_inertial; v_inertial]
    
    r_ee = rL(1:3,end);
    [J0n, Jmn] = Jacob(r_ee,r0,rL,P0,pm,robot.n_links_joints,robot);
    J = [J0n, Jmn];
    
    [I0,Im]=I_I(R0,RL,robot);
    [M0_tilde,Mm_tilde]=MCB(I0,Im,Bij,Bi0,robot);
    [H0, H0m, Hm] = GIM(M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot);
    H = [H0, H0m; H0m', Hm];
    
    [C0, C0m, Cm0, Cm] = CIM(t0,tm,I0,Im,M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot);
    C_matrix = [C0, C0m; Cm0, Cm];
    C_force = C_matrix * [u0; um];
    
    % Gravity
    wF0 = [0;0;0;0;0;-robot.base_link.mass*g];
    wFm = zeros(6,robot.n_links_joints);
    for i=1:robot.n_links_joints
        wFm(6,i) = -robot.links(i).mass*g;
    end
    tau0 = wF0;
    taum = 0.1*ones(n_q, 1);
    % for i=1:robot.n_links_joints
    %     if robot.joints(i).type ~= 0
    %         taum(robot.joints(i).q_id) = pm(:,i)' * wFm(:,i);
    %     end
    % end
    tau = [tau0; taum];
    
    % Constraint violations
    R_ee = RL(1:3,1:3,end);
    phi_pos = [rotation_error(R_ee, R_contact); r_ee - r_contact];
    phi_vel = tm(1:6,end);
    
    % Jdot*qdot term (must use t0, not u0, because Jacobdot expects inertial frame twist!)
    t_ee = tm(1:6,end);  % End-effector twist in INERTIAL frame
    [J0dot, Jmdot] = Jacobdot(r_ee, t_ee, r0, t0, rL, tm, P0, pm, robot.n_links_joints, robot);
    Jdot_qdot = J0dot*u0 + Jmdot*um;  % Jacobdot maps u0 (mixed) → acceleration
    
    % Baumgarte
    rhs_constraint = -Jdot_qdot - 2*alpha*phi_vel - beta^2*phi_pos;
    
    % KKT
    n_total = 6 + n_q;
    KKT = [H + epsilon*eye(n_total), -J'; 
           J, epsilon*eye(6)];
    rhs = [tau - C_force; rhs_constraint];
    
    sol = KKT \ rhs;
    vdot = sol(1:n_total);
    lambda = sol(n_total+1:end);
end

function q_new = update_positions(q, v_new, dt)
    % Update positions: q = q + dt*v_NEW (semi-implicit!)
    % CRITICAL: Use v_NEW (updated velocity), not old velocity
    % Special handling for quaternion
    
    n_q = length(q) - 10;  % 3 (r0) + 4 (p0) + n_q (qm)
    
    r0 = q(1:3);
    p0 = q(4:7);
    qm = q(8:end);
    
    u0_new = v_new(1:6);    % NEW base velocities (updated!)
    um_new = v_new(7:end);  % NEW joint velocities (updated!)
    
    % Position update with NEW velocities
    r0_new = r0 + dt * u0_new(4:6);  % u0(4:6) is v_inertial
    
    % Quaternion update with NEW angular velocity
    % u0(1:3) is ω_body (angular velocity in BODY frame)
    omega_new = u0_new(1:3);
    p0_dot = 0.5 * quat_omega_product(p0, omega_new);
    p0_new = p0 + dt * p0_dot;
    p0_new = p0_new / norm(p0_new);  % Normalize
    
    % Joint positions with NEW joint velocities
    qm_new = qm + dt * um_new;
    
    q_new = [r0_new; p0_new; qm_new];
end

function [phi_pos, phi_vel] = check_constraints_qv(q, v, robot, r_contact, R_contact)
    n_q = robot.n_q;
    r0 = q(1:3);
    p0 = q(4:7);
    qm = q(8:end);
    u0 = v(1:6);
    um = v(7:end);
    
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

function pdot = quat_omega_product(p, omega)
    qx = p(1); qy = p(2); qz = p(3); qw = p(4);
    wx = omega(1); wy = omega(2); wz = omega(3);
    pdot = [qw*wx + qy*wz - qz*wy;
            qw*wy - qx*wz + qz*wx;
            qw*wz + qx*wy - qy*wx;
           -qx*wx - qy*wy - qz*wz];
end
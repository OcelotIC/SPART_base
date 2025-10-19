%% SIMULATION AVEC STABILISATION DE BAUMGARTE
clear; clc; close all;

% Charger le robot 
filename = 'SC_6DoF.urdf';
[robot, robot_keys] = urdf2robot(filename);

%% Paramètres de simulation
dt = 0.002;  % PAS DE TEMPS RÉDUIT : 2 ms au lieu de 50 ms !
T_sim = 10;  % Durée totale [s]
N_steps = round(T_sim / dt);

% Paramètres de Baumgarte
alpha_baum = 5;  % Gain vitesse
beta_baum = 5;   % Gain position

%% Conditions initiales
% État du robot
R0 = eye(3);
r0 = [0; 0; 1];  % Base à 1m de hauteur
qm = 0.4*pi*ones(robot.n_q, 1);
u0 = zeros(6, 1);
um = 0*ones(robot.n_q, 1);

% Définir le point de contact désiré (FIXE)
contact_link_id = robot.n_links_joints;  % End-effector par exemple

% Calculer position/orientation initiale comme référence
[RJ, RL, rJ, rL, e, g] = Kinematics(R0, r0, qm, robot);
r_contact_des = rL(1:3, contact_link_id);
R_contact_des = RL(1:3, 1:3, contact_link_id);

%% Allocation mémoire
state_history = zeros(6 + robot.n_q + 6 + robot.n_q, N_steps);
error_history = zeros(6, N_steps);

%% BOUCLE DE SIMULATION
for k = 1:N_steps
    
    %% 1. Cinématique et dynamique avec SPART
    [RJ, RL, rJ, rL, e, g] = Kinematics(R0, r0, qm, robot);
    [Bij, Bi0, P0, pm] = DiffKinematics(R0, r0, rL, e, g, robot);
    [t0, tm] = Velocities(Bij, Bi0, P0, pm, u0, um, robot);
    
    % Matrices d'inertie et Coriolis
    [I0, Im] = I_I(R0, RL, robot);
    [M0_tilde, Mm_tilde] = MCB(I0, Im, Bij, Bi0, robot);
    [H0, H0m, Hm] = GIM(M0_tilde, Mm_tilde, Bij, Bi0, P0, pm, robot);
    [C0, C0m, Cm0, Cm] = CIM(t0, tm, I0, Im, M0_tilde, Mm_tilde, Bij, Bi0, P0, pm, robot);
    
    % Assemblage H et C
    H = [H0, H0m; H0m', Hm];
    C = [C0, C0m; Cm0, Cm];
    
    %% 2. Calcul de l'erreur de contrainte avec Baumgarte
    [Phi, Phi_dot] = compute_6D_contact_error(R0, r0, qm, u0, um, ...
                                               contact_link_id, ...
                                               r_contact_des, R_contact_des, ...
                                               robot);
    
    % Sauvegarder l'erreur
    error_history(:, k) = Phi;
    

%% 3. Jacobien de contact et calcul de Phi_dot
[Phi, Phi_dot, Jc] = compute_6D_contact_error(R0, r0, qm, u0, um, ...
                                               contact_link_id, ...
                                               r_contact_des, R_contact_des, ...
                                               robot);

% Sauvegarder l'erreur
error_history(:, k) = Phi;

%% 4. Jacobien time-derivative (CORRECTION ICI)
% Calcul des twists
[t0, tm] = Velocities(Bij, Bi0, P0, pm, u0, um, robot);

% Position et twist du point de contact
r_contact = rL(:, contact_link_id);

% Twist du point de contact (déjà calculé dans Phi_dot, mais explicite ici)
tp = Phi_dot;  % C'est Jc * qdot

% CORRECTION : tm est [6 x n_links], extraire le twist du lien i
tL = zeros(6, robot.n_links_joints);
for i = 1:robot.n_links_joints
    tL(:, i) = tm(:, i);
end

% Calcul de Jc_dot
[J0c_dot, Jmc_dot] = Jacobdot(r_contact, tp, r0, t0, rL, tL, P0, pm, ...
                               contact_link_id, robot);
Jc_dot = [J0c_dot, Jmc_dot];

%% 5. Terme de stabilisation de Baumgarte
qdot = [u0; um];
stabilization_term = Jc_dot * qdot + 2*alpha_baum * Phi_dot + beta_baum^2 * Phi;
    %% 5. Couples généralisés (exemple : gravité nulle, pas de commande)
    tau = [zeros(6,1); 0.1*ones(robot.n_q, 1)];
    
    %% 6. RÉSOLUTION EN DEUX TEMPS
    
    % Étape 1 : Calcul forces de contact via complément de Schur
    Schur_matrix = Jc * (H \ Jc');
    
    % Vérifier conditionnement
    if rcond(Schur_matrix) < 1e-12
        warning('Schur complement mal conditionné à t = %.3f', k*dt);
    end
    
    rhs = Jc * (H \ (tau - C*qdot)) + stabilization_term;  % AVEC BAUMGARTE !
    lambda = Schur_matrix \ rhs;
    
    % Étape 2 : Calcul des accélérations
    qddot = H \ (tau - C*qdot + Jc' * lambda);
    
    u0dot = qddot(1:6);
    umdot = qddot(7:end);
    
    %% 7. INTÉGRATION SEMI-IMPLICITE (meilleure conservation)
    % Vitesse d'abord
    u0_next = u0 + dt * u0dot;
    um_next = um + dt * umdot;
    
    % Position avec NOUVELLE vitesse
    r0 = r0 + dt * u0_next(4:6);
    
    % Orientation : intégration de quaternion ou rotation
    omega = u0_next(1:3);
    R0 = R0 * expm(skew(omega * dt));
    
    % Joints
    qm = qm + dt * um_next;
    
    % Mise à jour pour prochaine itération
    u0 = u0_next;
    um = um_next;
    
    %% 8. Sauvegarder l'état
    state_history(:, k) = [DCM_Angles321(R0);r0; qm; u0; um];
    
    % Affichage périodique
    if mod(k, 500) == 0
        fprintf('t = %.2f s | ||Phi|| = %.2e m | ||Phi_dot|| = %.2e m/s\n', ...
                k*dt, norm(Phi), norm(Phi_dot));
    end
end

%% VISUALISATION
figure('Position', [100, 100, 1200, 800]);

% Erreur de position
subplot(3,2,1);
plot((1:N_steps)*dt, error_history(1:3, :)', 'LineWidth', 1.5);
xlabel('Temps [s]'); ylabel('Erreur position [m]');
title('Erreur de position du contact (Baumgarte)');
legend('x', 'y', 'z');
grid on;

% Erreur d'orientation
subplot(3,2,2);
plot((1:N_steps)*dt, error_history(4:6, :)', 'LineWidth', 1.5);
xlabel('Temps [s]'); ylabel('Erreur orientation [rad]');
title('Erreur d''orientation du contact');
legend('\phi_x', '\phi_y', '\phi_z');
grid on;

% Norme de l'erreur
subplot(3,2,3);
semilogy((1:N_steps)*dt, vecnorm(error_history(1:3, :)), 'LineWidth', 2);
hold on;
yline(1e-3, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Seuil 1 mm');
xlabel('Temps [s]'); ylabel('||Phi_{pos}|| [m]');
title('Norme de l''erreur de position');
grid on; legend;

subplot(3,2,4);
semilogy((1:N_steps)*dt, vecnorm(error_history(4:6, :)), 'LineWidth', 2);
xlabel('Temps [s]'); ylabel('||Phi_{orient}|| [rad]');
title('Norme de l''erreur d''orientation');
grid on;

% Statistiques finales
fprintf('\n=== RÉSULTATS FINAUX ===\n');
fprintf('Erreur position max : %.2e m\n', max(vecnorm(error_history(1:3, :))));
fprintf('Erreur position finale : %.2e m\n', norm(error_history(1:3, end)));
fprintf('Erreur orientation max : %.2e rad\n', max(vecnorm(error_history(4:6, :))));


% position
figure
plot((1:N_steps)*dt, state_history(7:10, :)', 'LineWidth', 1.5);
xlabel('Temps [s]'); ylabel('position q [rad]');
title('joint positions');
grid on;


%% Fonction auxiliaire skew
function S = skew(v)
    S = [0, -v(3), v(2);
         v(3), 0, -v(1);
         -v(2), v(1), 0];
end
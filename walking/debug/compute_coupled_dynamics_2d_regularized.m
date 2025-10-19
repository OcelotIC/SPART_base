function dstate = compute_coupled_dynamics_2d_regularized(t, state_vec, walker, satellite, params)
% COMPUTE_COUPLED_DYNAMICS_2D_REGULARIZED - Dynamique couplée avec régularisation
%
% Version stabilisée avec:
%   - Régularisation Tikhonov de la matrice de masse augmentée
%   - Détection géométrique des contacts (pas de contraintes fantômes)
%   - Gestion robuste des singularités
%   - Projection 2D propre et conservatrice
%
% Inputs:
%   t         : Temps [s]
%   state_vec : Vecteur d'état [n×1]
%   walker    : Structure robot walker
%   satellite : Structure robot satellite
%   params    : Paramètres de régularisation
%       .lambda_tikhonov : Régularisation (défaut: 1e-6)
%       .contact_threshold : Distance contact (défaut: 0.05 m)
%       .cond_max : Cond. max acceptable (défaut: 1e8)
%
% Output:
%   dstate : Dérivée temporelle de state_vec

%% PARAMÈTRES DE RÉGULARISATION
if nargin < 5
    params = struct();
end
if ~isfield(params, 'lambda_tikhonov')
    params.lambda_tikhonov = 1e-6; % Régularisation Tikhonov
end
if ~isfield(params, 'contact_threshold')
    params.contact_threshold = 0.05; % 5 cm
end
if ~isfield(params, 'cond_max')
    params.cond_max = 1e8;
end

%% 1. UNPACKING DE L'ÉTAT
n_walker = 6 + 6; % 6 DoF base + 6 joints
n_sat = 6 + 1;    % 6 DoF base + 1 RW

idx = 1;

% Walker
r0_w = state_vec(idx:idx+2); idx = idx + 3;
quat_w = state_vec(idx:idx+3); idx = idx + 4;
qm_w = state_vec(idx:idx+5); idx = idx + 6;
u0_w = state_vec(idx:idx+5); idx = idx + 6;
um_w = state_vec(idx:idx+5); idx = idx + 6;

% Satellite
r0_s = state_vec(idx:idx+2); idx = idx + 3;
quat_s = state_vec(idx:idx+3); idx = idx + 4;
u0_s = state_vec(idx:idx+5); idx = idx + 6;
h_RW = state_vec(idx); 

%% 2. NORMALISATION DES QUATERNIONS
quat_w = quat_w / norm(quat_w);
quat_s = quat_s / norm(quat_s);

%% 3. PROJECTION 2D (x, z, θ_y)
P_2d = zeros(3, 6);
P_2d(1, 1) = 1;  % Translation X
P_2d(2, 3) = 1;  % Translation Z
P_2d(3, 5) = 1;  % Rotation Y

% États 2D projetés
r0_w_2d = P_2d * [r0_w; 0; 0; 0];  % [x; z; 0]
r0_s_2d = P_2d * [r0_s; 0; 0; 0];

u0_w_2d = P_2d * u0_w;  % [vx; vz; ωy]
u0_s_2d = P_2d * u0_s;

%% 4. DÉTECTION GÉOMÉTRIQUE DES CONTACTS
[R0_w, ~] = quat2DCM(quat_w);
[R0_s, ~] = quat2DCM(quat_s);

% Positions des end-effectors (via FK rapide)
[r_ee1, r_ee2] = compute_ee_positions_fast(r0_w, R0_w, qm_w);

% Positions des docks
dock1_local = [1.5; 0; 0];
dock2_local = [-1.5; 0; 0];
r_dock1 = r0_s + R0_s * dock1_local;
r_dock2 = r0_s + R0_s * dock2_local;

% Distances
d1 = norm(r_ee1 - r_dock1);
d2 = norm(r_ee2 - r_dock2);

% État des contacts
contact1_active = (d1 < params.contact_threshold);
contact2_active = (d2 < params.contact_threshold);
n_contacts_active = contact1_active + contact2_active;

%% 5. MATRICE DE MASSE (GIM via SPART)
[~, ~, ~, M_walker] = FD(walker, r0_w, quat_w, qm_w, zeros(6,1), zeros(6,1), zeros(6,1));
[~, ~, ~, M_sat] = FD(satellite, r0_s, quat_s, 0, zeros(6,1), 0, zeros(6,1));

% Projection 2D des matrices de masse
M_w_base_2d = P_2d * M_walker(1:6, 1:6) * P_2d';  % 3×3
M_w_joints = M_walker(7:12, 7:12);                % 6×6

M_s_base_2d = P_2d * M_sat(1:6, 1:6) * P_2d';    % 3×3
I_RW = satellite.I_RW;                            % Scalaire

%% 6. JACOBIENNE DE CONTACT (si contacts actifs)
if n_contacts_active > 0
    n_constraints = 2 * n_contacts_active; % 2 contraintes (x,z) par contact
    
    % Taille totale : [3_w + 6_w + 3_s + 1_RW] = 13
    Jc = zeros(n_constraints, 13);
    
    constraint_idx = 1;
    
    % Contact 1
    if contact1_active
        % Jacobienne end-effector 1 (simplifié - à remplacer par GIM_Jacobian)
        J_ee1_full = compute_jacobian_ee(walker, qm_w, 1); % 6×12 (base + joints)
        J_ee1_2d = [P_2d * J_ee1_full(1:3, 1:6), J_ee1_full(1:3, 7:12)]; % [3_base 6_joints]
        
        % Contraintes X et Z seulement
        Jc(constraint_idx, 1:9) = J_ee1_2d([1,3], :); % Lignes X et Z
        constraint_idx = constraint_idx + 2;
    end
    
    % Contact 2
    if contact2_active
        J_ee2_full = compute_jacobian_ee(walker, qm_w, 2);
        J_ee2_2d = [P_2d * J_ee2_full(1:3, 1:6), J_ee2_full(1:3, 7:12)];
        
        Jc(constraint_idx, 1:9) = J_ee2_2d([1,3], :);
        constraint_idx = constraint_idx + 2;
    end
    
    % Satellite contribue aussi au mouvement relatif
    Jc(:, 10:12) = -Jc(:, 1:3); % Position satellite opposée
    
else
    % Pas de contacts actifs
    n_constraints = 0;
    Jc = [];
end

%% 7. ASSEMBLAGE MATRICE DE MASSE AUGMENTÉE
M_sys = blkdiag(M_w_base_2d, M_w_joints, M_s_base_2d, I_RW);

if n_constraints > 0
    % Avec contraintes
    M_aug = [M_sys, -Jc';
             Jc, zeros(n_constraints, n_constraints)];
    
    % RÉGULARISATION TIKHONOV
    reg_matrix = params.lambda_tikhonov * eye(size(M_aug));
    reg_matrix(1:13, 1:13) = 0; % Ne régulariser que le bloc des multiplicateurs
    M_aug_reg = M_aug + reg_matrix;
    
    % Vérification conditionnement
    cond_M = cond(M_aug_reg);
    if cond_M > params.cond_max
        warning('t=%.3f: Conditionnement élevé (%.2e) - augmentation régularisation', t, cond_M);
        % Régularisation adaptative
        M_aug_reg = M_aug + (10 * params.lambda_tikhonov) * eye(size(M_aug));
    end
    
else
    % Sans contraintes (phase swing)
    M_aug_reg = M_sys;
end

%% 8. FORCES GÉNÉRALISÉES
% Gravité (2D: seulement en Z)
g_2d = [0; -9.81; 0];

% Walker
f_grav_w_base = M_w_base_2d * g_2d;
f_grav_w_joints = zeros(6, 1); % Pas de gravité sur les joints

% Satellite (flottant - pas de gravité orbitale)
f_grav_s = zeros(3, 1);

% Couple de contrôle RW (pour l'instant nul)
tau_RW = 0;

% Vecteur des forces généralisées
Q_sys = [f_grav_w_base; f_grav_w_joints; f_grav_s; tau_RW];

if n_constraints > 0
    Q_aug = [Q_sys; zeros(n_constraints, 1)]; % Pas de forces aux contacts (holonomiques)
else
    Q_aug = Q_sys;
end

%% 9. RÉSOLUTION DU SYSTÈME
% M_aug * [accel; lambda] = Q_aug
try
    sol = M_aug_reg \ Q_aug;
catch ME
    warning('t=%.3f: Échec résolution linéaire - %s', t, ME.message);
    % Retour sécurisé : vitesses constantes
    dstate = zeros(size(state_vec));
    return;
end

% Extraction des accélérations
accel_w_base = sol(1:3);
accel_w_joints = sol(4:9);
accel_s_base = sol(10:12);
accel_h_RW = sol(13);

% Forces de contact (pour monitoring)
if n_constraints > 0
    lambda = sol(14:end);
    % Stockage pour post-traitement (optionnel)
else
    lambda = [];
end

%% 10. DÉTECTION D'INSTABILITÉS
if any(~isfinite(sol))
    warning('t=%.3f: NaN/Inf détecté dans les accélérations!', t);
    dstate = zeros(size(state_vec));
    return;
end

if norm(accel_w_joints) > 1000 % Accélérations irréalistes
    warning('t=%.3f: Accélérations articulaires explosives (%.2e rad/s²)', t, norm(accel_w_joints));
end

%% 11. ASSEMBLAGE DÉRIVÉE D'ÉTAT
dstate = zeros(size(state_vec));

idx = 1;

% Walker
dstate(idx:idx+2) = u0_w_2d; idx = idx + 3;
dstate(idx:idx+3) = 0.5 * G(quat_w) * u0_w; idx = idx + 4; % Dérivée quaternion
dstate(idx:idx+5) = um_w; idx = idx + 6;
dstate(idx:idx+5) = [accel_w_base; zeros(3,1)]; idx = idx + 6; % Expansion 2D→6D
dstate(idx:idx+5) = accel_w_joints; idx = idx + 6;

% Satellite
dstate(idx:idx+2) = u0_s_2d; idx = idx + 3;
dstate(idx:idx+3) = 0.5 * G(quat_s) * u0_s; idx = idx + 4;
dstate(idx:idx+5) = [accel_s_base; zeros(3,1)]; idx = idx + 6;
dstate(idx) = accel_h_RW;

end

%% ============================================================================
%% FONCTIONS AUXILIAIRES
%% ============================================================================

function [r_ee1, r_ee2] = compute_ee_positions_fast(r0, R0, qm)
% FK rapide pour end-effectors (sans SPART)
L = 0.4;

% Bras 1
r_attach1 = r0 + R0 * [0.3; 0; 0];
q1 = qm(1:3);
x1 = L * (cos(q1(1)) + cos(q1(1)+q1(2)) + cos(q1(1)+q1(2)+q1(3)));
z1 = L * (sin(q1(1)) + sin(q1(1)+q1(2)) + sin(q1(1)+q1(2)+q1(3)));
r_ee1 = r_attach1 + R0 * [x1; 0; z1];

% Bras 2
r_attach2 = r0 + R0 * [-0.3; 0; 0];
q2 = qm(4:6);
x2 = L * (cos(q2(1)) + cos(q2(1)+q2(2)) + cos(q2(1)+q2(2)+q2(3)));
z2 = L * (sin(q2(1)) + sin(q2(1)+q2(2)) + sin(q2(1)+q2(2)+q2(3)));
r_ee2 = r_attach2 + R0 * [x2; 0; z2];
end

function J = compute_jacobian_ee(robot, qm, arm_id)
% Jacobienne end-effector (simplifié - à remplacer par GIM_Jacobian SPART)
% Pour l'instant : approximation numérique

L = 0.4;

if arm_id == 1
    q = qm(1:3);
    offset = [0.3; 0; 0];
else
    q = qm(4:6);
    offset = [-0.3; 0; 0];
end

% Jacobienne analytique 3R planaire
J_arm = zeros(6, 12);

% Contribution de la base (translation + rotation)
J_arm(1:3, 1:3) = eye(3);
J_arm(1:3, 4:6) = -skew(offset); % Effet rotation base

% Contribution des joints
for i = 1:3
    % Colonne pour joint i
    dtheta = sum(q(1:i));
    dL = L * i;
    
    J_arm(1, 6 + (arm_id-1)*3 + i) = -dL * sin(dtheta);
    J_arm(3, 6 + (arm_id-1)*3 + i) = dL * cos(dtheta);
end

J = J_arm;
end

function S = skew(v)
% Matrice antisymétrique
S = [0 -v(3) v(2);
     v(3) 0 -v(1);
     -v(2) v(1) 0];
end

function G_mat = G(q)
% Matrice G pour dérivée de quaternion
% q_dot = 0.5 * G(q) * omega
epsilon = q(1:3);
eta = q(4);
G_mat = [-skew(epsilon) + eta*eye(3); -epsilon'];
end
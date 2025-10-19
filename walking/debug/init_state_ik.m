function state = initialize_state_ik(walker, satellite, config)
% INITIALIZE_STATE_IK - Crée un état initial géométriquement cohérent
%
% Calcule les angles articulaires pour que les end-effectors touchent
% EXACTEMENT les docks du satellite (via IK)
%
% Inputs:
%   walker    - Structure robot walker
%   satellite - Structure robot satellite  
%   config    - Structure avec:
%       .r0_sat       : Position satellite [3×1]
%       .theta_sat_2d : Rotation satellite autour Y (rad)
%       .z_walker     : Hauteur du torse walker (m)
%       .H_initial    : Momentum initial [3×1] (optionnel)
%
% Output:
%   state - État initial compatible avec compute_coupled_dynamics_2d

fprintf('\n╔════════════════════════════════════════════╗\n');
fprintf('║   INITIALISATION PAR CINÉMATIQUE INVERSE   ║\n');
fprintf('╚════════════════════════════════════════════╝\n\n');

%% 1. CONFIGURATION SATELLITE
state.r0_sat = config.r0_sat;

% Quaternion depuis rotation 2D (rotation autour Y)
theta = config.theta_sat_2d;
state.quat_sat = [0; sin(theta/2); 0; cos(theta/2)];

state.u0_sat = zeros(6, 1);
state.h_RW = 0;

[R0_sat] = quat_DCM(state.quat_sat);

fprintf('Satellite:\n');
fprintf('  Position: [%.2f, %.2f, %.2f] m\n', state.r0_sat);
fprintf('  Rotation: %.1f°\n\n', rad2deg(theta));

%% 2. POSITIONS DES DOCKS (CIBLES IK)
% Suppose docks à ±1.5m en X du centre satellite (à adapter selon votre modèle)
dock1_local = [1.5; 0; 0];   % Dock droit
dock2_local = [-1.5; 0; 0];  % Dock gauche

r_dock1_world = state.r0_sat + R0_sat * dock1_local;
r_dock2_world = state.r0_sat + R0_sat * dock2_local;

fprintf('Docks (repère inertiel):\n');
fprintf('  Dock 1: [%.2f, %.2f, %.2f] m\n', r_dock1_world);
fprintf('  Dock 2: [%.2f, %.2f, %.2f] m\n\n', r_dock2_world);

%% 3. CONFIGURATION WALKER - BASE
% Position du torse: entre les deux docks, à la hauteur spécifiée
state.r0_walker = [(r_dock1_world(1) + r_dock2_world(1))/2;
                   0;
                   config.z_walker];

% Orientation: alignée avec satellite en 2D
state.quat_walker = state.quat_sat;

state.u0_walker = zeros(6, 1);
state.um_walker = zeros(6, 1);

[R0_walker] = quat_DCM(state.quat_walker);

fprintf('Walker torse:\n');
fprintf('  Position: [%.2f, %.2f, %.2f] m\n', state.r0_walker);
fprintf('  Orientation: %.1f°\n\n', rad2deg(theta));

%% 4. INVERSE KINEMATICS - BRAS 1 (DROIT)
% Position de l'attache du bras 1 sur le torse
r_attach1_local = [0.3; 0; 0]; % 30 cm à droite du centre torse
r_attach1 = state.r0_walker + R0_walker * r_attach1_local;

% Vecteur attache → dock dans le repère du bras
r_target1_arm_frame = R0_walker' * (r_dock1_world - r_attach1);

% IK planaire 3R (plan XZ local)
L = 0.4; % Longueur des liens
q1 = solve_ik_3r_planar(r_target1_arm_frame([1,3]), L);

if isempty(q1)
    error('IK infaisable pour bras 1 (cible hors workspace)');
end

fprintf('Bras 1 (IK):\n');
fprintf('  Angles: [%.1f°, %.1f°, %.1f°]\n', rad2deg(q1));

% Vérification FK
r_ee1_check = forward_kinematics_3r(q1, L);
r_ee1_world = r_attach1 + R0_walker * [r_ee1_check(1); 0; r_ee1_check(2)];
error1 = norm(r_ee1_world - r_dock1_world);
fprintf('  Erreur position: %.4f m\n', error1);

if error1 > 1e-3
    warning('IK bras 1: erreur > 1mm');
end

%% 5. INVERSE KINEMATICS - BRAS 2 (GAUCHE)
r_attach2_local = [-0.3; 0; 0];
r_attach2 = state.r0_walker + R0_walker * r_attach2_local;

r_target2_arm_frame = R0_walker' * (r_dock2_world - r_attach2);

q2 = solve_ik_3r_planar(r_target2_arm_frame([1,3]), L);

if isempty(q2)
    error('IK infaisable pour bras 2 (cible hors workspace)');
end

fprintf('\nBras 2 (IK):\n');
fprintf('  Angles: [%.1f°, %.1f°, %.1f°]\n', rad2deg(q2));

r_ee2_check = forward_kinematics_3r(q2, L);
r_ee2_world = r_attach2 + R0_walker * [r_ee2_check(1); 0; r_ee2_check(2)];
error2 = norm(r_ee2_world - r_dock2_world);
fprintf('  Erreur position: %.4f m\n', error2);

if error2 > 1e-3
    warning('IK bras 2: erreur > 1mm');
end

%% 6. ASSEMBLAGE VECTEUR ARTICULAIRE
state.qm_walker = [q1; q2];

%% 7. MOMENTUM INITIAL (SI SPÉCIFIÉ)
if isfield(config, 'H_initial')
    state.h_total_initial = config.H_initial;
else
    state.h_total_initial = [0; 0; 0]; % Système au repos
end

%% 8. RÉSUMÉ
fprintf('\n═══════════════════════════════════════════\n');
fprintf('✅ Initialisation terminée:\n');
fprintf('   - End-effectors sur docks (err < %.1e m)\n', max(error1, error2));
fprintf('   - Configuration articulaire calculée par IK\n');
fprintf('   - Système au repos (vitesses nulles)\n');
fprintf('═══════════════════════════════════════════\n\n');

fprintf('⚠️  Lancez maintenant: diagnose_initial_state(state, walker, satellite)\n\n');

end

%% ============================================================================
%% FONCTIONS AUXILIAIRES
%% ============================================================================

function q = solve_ik_3r_planar(target, L)
% IK pour chaîne 3R planaire (3 liens égaux dans plan XZ)
%
% Inputs:
%   target : Position cible [x; z] dans repère du bras
%   L      : Longueur de chaque lien
%
% Output:
%   q : Angles [q1; q2; q3] ou [] si infaisable
%
% Méthode: Géométrique avec redondance (plusieurs solutions possibles)

x = target(1);
z = target(2);
d = sqrt(x^2 + z^2);

% Vérification workspace
L_total = 3 * L;
if d > L_total || d < 0.1 % Minimum pour éviter singularités
    q = [];
    return;
end

% Solution analytique simplifiée (coude "up")
% On cherche une config avec q2 négatif (coude vers le bas)

% Angle au poignet (lien 3)
psi = atan2(z, x);

% Distance projetée
d_proj = d - L; % Compenser le dernier lien

% Loi des cosinus pour 2R
c2 = (d_proj^2 - 2*L^2) / (2*L^2);
c2 = max(-1, min(1, c2)); % Clamp pour stabilité numérique

q2 = -acos(c2); % Coude vers le bas

% Angle du premier lien
alpha = atan2(L*sin(q2), L + L*cos(q2));
q1 = psi - alpha;

% Angle du troisième lien (pour atteindre la cible)
q3 = -(q1 + q2); % Force la colinéarité finale

q = [q1; q2; q3];

% Vérification (sécurité)
r_check = forward_kinematics_3r(q, L);
err = norm([x; z] - r_check);
if err > 1e-2
    warning('IK: erreur résiduelle = %.4f m', err);
end

end

function r = forward_kinematics_3r(q, L)
% Forward kinematics 3R planaire
%
% Inputs:
%   q : Angles [q1; q2; q3]
%   L : Longueur des liens
%
% Output:
%   r : Position end-effector [x; z]

theta1 = q(1);
theta2 = q(1) + q(2);
theta3 = q(1) + q(2) + q(3);

x = L * (cos(theta1) + cos(theta2) + cos(theta3));
z = L * (sin(theta1) + sin(theta2) + sin(theta3));

r = [x; z];
end
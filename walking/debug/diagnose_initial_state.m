function diagnostics = diagnose_initial_state(state, walker, satellite)
% DIAGNOSE_INITIAL_STATE - Vérifie la cohérence de l'état initial
%
% Identifie les problèmes numériques AVANT de lancer la simulation
%
% Usage:
%   diag = diagnose_initial_state(state, walker, satellite);
%   if ~diag.is_valid
%       error('État initial invalide : %s', diag.error_msg);
%   end

fprintf('\n╔════════════════════════════════════════════╗\n');
fprintf('║   DIAGNOSTIC DE L''ÉTAT INITIAL 2D          ║\n');
fprintf('╚════════════════════════════════════════════╝\n\n');

diagnostics = struct();
diagnostics.is_valid = true;
diagnostics.error_msg = '';
diagnostics.warnings = {};

%% 1. VÉRIFICATION GÉOMÉTRIQUE DES CONTACTS
fprintf('1. Géométrie des contacts\n');
fprintf('   -------------------------\n');

% Position des end-effectors dans le repère inertiel
[R0_walker, ~] = quat2DCM(state.quat_walker);
r0_walker = state.r0_walker;

% Forward kinematics pour trouver les positions des EE
% (Simplifié - suppose que les EE sont les derniers links)
% NOTE: Adapter selon votre modèle exact

% Bras 1 (links 1-3)
q1 = state.qm_walker(1:3);
r_ee1 = compute_ee_position(r0_walker, R0_walker, q1, walker, 'arm1');

% Bras 2 (links 4-6)
q2 = state.qm_walker(4:6);
r_ee2 = compute_ee_position(r0_walker, R0_walker, q2, walker, 'arm2');

% Position des docks sur le satellite
[R0_sat, ~] = quat2DCM(state.quat_sat);
r0_sat = state.r0_sat;

% Docks positions (à adapter selon votre modèle satellite)
% Suppose docks à ±1.5m en x du centre satellite
r_dock1 = r0_sat + R0_sat * [1.5; 0; 0];
r_dock2 = r0_sat + R0_sat * [-1.5; 0; 0];

% Distances EE-Dock
d1 = norm(r_ee1 - r_dock1);
d2 = norm(r_ee2 - r_dock2);

fprintf('   Distance EE1 → Dock1: %.4f m\n', d1);
fprintf('   Distance EE2 → Dock2: %.4f m\n', d2);

tolerance_contact = 0.05; % 5 cm
if d1 > tolerance_contact
    diagnostics.warnings{end+1} = sprintf('EE1 loin du dock (%.3f m)', d1);
    fprintf('   ⚠️  EE1 n''est pas en contact!\n');
end
if d2 > tolerance_contact
    diagnostics.warnings{end+1} = sprintf('EE2 loin du dock (%.3f m)', d2);
    fprintf('   ⚠️  EE2 n''est pas en contact!\n');
end

fprintf('\n');

%% 2. CALCUL DE LA MATRICE DE MASSE AUGMENTÉE
fprintf('2. Conditionnement de la matrice de masse\n');
fprintf('   ---------------------------------------\n');

% Matrice de masse du système couplé
n_walker = 6 + 6; % 6 base + 6 joints
n_sat = 6 + 1;    % 6 base + 1 RW

% Masses individuelles (simplifié)
M_walker_base = walker.base_link.inertia;
M_sat_base = satellite.base_link.inertia;

% Projection 2D (x, z, θ_y)
P_2d = zeros(3, 6);
P_2d(1, 1) = 1;  % x
P_2d(2, 3) = 1;  % z
P_2d(3, 5) = 1;  % θ_y

% Matrice de masse simplifiée (sans GIM complet pour diagnostic rapide)
M_w_2d = P_2d * M_walker_base(4:6, 4:6) * P_2d';
M_s_2d = P_2d * M_sat_base(4:6, 4:6) * P_2d';

% Jacobienne de contact (2D simplifié)
% Contact en 2D: 2 contraintes par point (x, z)
n_contacts = 2; % 2 end-effectors
n_constraints = 2 * n_contacts; % 4 contraintes (x,z pour chaque EE)

% Jacobienne approximative (à remplacer par compute_contact_jacobian_2d)
Jc = zeros(n_constraints, 3 + 6 + 3 + 1); % [walker_base walker_joints sat_base sat_RW]

% Matrice de masse augmentée
M_aug = [M_w_2d, zeros(3, 6), zeros(3, 3), zeros(3, 1);
         zeros(6, 3), eye(6)*5, zeros(6, 3), zeros(6, 1);  % Inertie articulaire ~5 kg·m²
         zeros(3, 3), zeros(3, 6), M_s_2d, zeros(3, 1);
         zeros(1, 3), zeros(1, 6), zeros(1, 3), satellite.I_RW];

M_aug = [M_aug, -Jc';
         Jc, zeros(n_constraints, n_constraints)];

% Conditionnement
cond_M = cond(M_aug);
fprintf('   Conditionnement M_aug: %.2e\n', cond_M);

if cond_M > 1e10
    diagnostics.is_valid = false;
    diagnostics.error_msg = sprintf('Matrice mal conditionnée (cond=%.2e)', cond_M);
    fprintf('   ❌ CRITIQUE: Matrice singulière!\n');
elseif cond_M > 1e6
    diagnostics.warnings{end+1} = sprintf('Conditionnement élevé (%.2e)', cond_M);
    fprintf('   ⚠️  Conditionnement élevé (risque instabilité)\n');
else
    fprintf('   ✅ Conditionnement acceptable\n');
end

% Valeurs propres
eigvals = eig(M_aug(1:end-n_constraints, 1:end-n_constraints));
min_eig = min(abs(eigvals));
fprintf('   Plus petite valeur propre: %.2e\n', min_eig);

if min_eig < 1e-8
    diagnostics.warnings{end+1} = 'Valeur propre quasi-nulle détectée';
    fprintf('   ⚠️  Valeur propre proche de zéro!\n');
end

fprintf('\n');

%% 3. CONSERVATION DU MOMENTUM
fprintf('3. Momentum initial\n');
fprintf('   ----------------\n');

% Vitesses nulles → H devrait être nul
H_walker = M_walker_base * [state.u0_walker(1:3); state.u0_walker(4:6)] + ...
           zeros(6, 1); % Contribution articulaire (vitesses nulles)

H_sat = M_sat_base * [state.u0_sat(1:3); state.u0_sat(4:6)] + ...
        [0; 0; 0; 0; 0; state.h_RW];

H_total_calc = H_walker + H_sat;
H_total_2d = P_2d * H_total_calc(4:6); % Projection en 2D

fprintf('   H_total (2D): [%.2e, %.2e, %.2e]\n', H_total_2d(1), H_total_2d(2), H_total_2d(3));
fprintf('   H_RW initial: %.2e N·m·s\n', state.h_RW);

if norm(H_total_2d) > 1e-6
    diagnostics.warnings{end+1} = 'Momentum initial non nul';
    fprintf('   ⚠️  Momentum initial devrait être nul!\n');
else
    fprintf('   ✅ Momentum conservé à t=0\n');
end

fprintf('\n');

%% 4. VÉRIFICATION DES LIMITES PHYSIQUES
fprintf('4. Limites physiques\n');
fprintf('   ------------------\n');

% Angles articulaires
qm_deg = rad2deg(state.qm_walker);
fprintf('   Angles bras 1: [%.1f°, %.1f°, %.1f°]\n', qm_deg(1:3));
fprintf('   Angles bras 2: [%.1f°, %.1f°, %.1f°]\n', qm_deg(4:6));

% Limites typiques: ±180°
if any(abs(state.qm_walker) > pi)
    diagnostics.warnings{end+1} = 'Angles articulaires hors limites';
    fprintf('   ⚠️  Angles > 180°!\n');
else
    fprintf('   ✅ Angles dans limites\n');
end

% Position walker (devrait être au-dessus du satellite)
if state.r0_walker(3) < state.r0_sat(3)
    diagnostics.warnings{end+1} = 'Walker sous le satellite';
    fprintf('   ⚠️  Walker sous le satellite!\n');
end

fprintf('\n');

%% 5. RÉSUMÉ
fprintf('═══════════════════════════════════════════\n');
fprintf('RÉSUMÉ:\n');
if diagnostics.is_valid
    fprintf('✅ État initial VALIDE pour simulation\n');
else
    fprintf('❌ État initial INVALIDE: %s\n', diagnostics.error_msg);
end

if ~isempty(diagnostics.warnings)
    fprintf('\n⚠️  Avertissements (%d):\n', length(diagnostics.warnings));
    for i = 1:length(diagnostics.warnings)
        fprintf('   %d. %s\n', i, diagnostics.warnings{i});
    end
end
fprintf('═══════════════════════════════════════════\n\n');

end

%% HELPER FUNCTIONS
function r_ee = compute_ee_position(r0_base, R0_base, q_arm, robot, arm_name)
% Calcule la position de l'end-effector par forward kinematics
% Simplifié pour 3 DoF planaires

L = 0.4; % Longueur des liens (à adapter)

% Position locale de l'attache du bras sur le torse
if strcmp(arm_name, 'arm1')
    r_attach_local = [0.3; 0; 0]; % Côté droit
else
    r_attach_local = [-0.3; 0; 0]; % Côté gauche
end

r_attach = r0_base + R0_base * r_attach_local;

% FK planaire (chaîne série 3R)
theta1 = q_arm(1);
theta2 = q_arm(1) + q_arm(2);
theta3 = q_arm(1) + q_arm(2) + q_arm(3);

% Position EE dans le plan XZ
x_ee_local = L * (cos(theta1) + cos(theta2) + cos(theta3));
z_ee_local = L * (sin(theta1) + sin(theta2) + sin(theta3));

r_ee = r_attach + R0_base * [x_ee_local; 0; z_ee_local];
end
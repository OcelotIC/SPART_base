function [robot, robot_keys] = create_walker_2d()
% CREATE_WALKER_2D - Creates a 2D walking robot from URDF
%
% Structure: Torso (floating base) + 2 arms × 3 links each
% Total: 9 DoF (3 base + 6 joints)
%
% Output:
%   robot      - SPART robot structure
%   robot_keys - Maps for link/joint IDs
%
% Uses SPART's urdf2robot() function (PROPER METHOD for branching structures)

%% URDF file path
% Assumes walker_2d.urdf is in the same directory or on path
urdf_file = 'walker_2d.urdf';

if ~exist(urdf_file, 'file')
    error(['URDF file not found: ', urdf_file, '\n', ...
           'Make sure walker_2d.urdf is in the current directory or MATLAB path.']);
end

%% Load robot from URDF using SPART
[robot, robot_keys] = urdf2robot(urdf_file);

%% Display summary
fprintf('✓ Walker 2D created from URDF\n');
fprintf('  Total mass: %.1f kg\n', robot.base_link.mass + sum([robot.links.mass]));
fprintf('  DoF: %d (6 joints)\n', robot.n_q);
fprintf('  Links: %d\n', robot.n_links_joints);

% Display connectivity
fprintf('  Connectivity (from SPART ConnectivityMap):\n');
fprintf('    child_base links: ');
child_base_indices = find(robot.con.child_base);
fprintf('%d ', child_base_indices);
fprintf('\n');

end
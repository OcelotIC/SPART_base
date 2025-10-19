function [robot_sat, satellite_params] = create_satellite_2d()
% CREATE_SATELLITE_2D - Creates a 2D satellite platform from URDF
%
% Structure: Rigid body (floating base) + 1 reaction wheel
% DoF: 3 (x, z, θ) + 1 RW momentum h_RW
%
% Output:
%   robot_sat        - SPART robot structure (from URDF)
%   satellite_params - Additional parameters (RW, contact points)
%
% Uses SPART's urdf2robot() function

%% URDF file path
urdf_file = 'satellite_2d.urdf';

if ~exist(urdf_file, 'file')
    error(['URDF file not found: ', urdf_file, '\n', ...
           'Make sure satellite_2d.urdf is in the current directory or MATLAB path.']);
end

%% Load satellite from URDF using SPART
[robot_sat, ~] = urdf2robot(urdf_file);

%% Satellite geometry parameters (for contact points and visualization)
L_sat = 5.0;   % Length [m]
W_sat = 2.0;   % Width [m]
H_sat = 0.3;   % Height [m]

%% Reaction Wheel Parameters
satellite_params.RW = struct();
satellite_params.RW.h_max = 30;      % Maximum momentum storage [N·m·s]
satellite_params.RW.tau_max = 3;     % Maximum torque [N·m]
satellite_params.RW.I_wheel = 0.5;   % Wheel inertia [kg·m²]
satellite_params.RW.axis = [0; 0; 1]; % Rotation around Z (2D)

%% Contact/Docking Points (local coordinates)
% Two docking interfaces on the satellite surface
% Positioned symmetrically at ±1.5m from center

satellite_params.contact_points = struct();

% Contact 1 (Right side, +X direction)
satellite_params.contact_points(1).position_local = [1.5; 0; 0.2];  % [x; y; z] in satellite frame
satellite_params.contact_points(1).orientation_local = 0;  % θ relative to satellite

% Contact 2 (Left side, -X direction)
satellite_params.contact_points(2).position_local = [-1.5; 0; 0.2];
satellite_params.contact_points(2).orientation_local = 0;

%% Geometric Properties (for visualization)
satellite_params.geometry = struct();
satellite_params.geometry.length = L_sat;
satellite_params.geometry.width = W_sat;
satellite_params.geometry.height = H_sat;

% Vertices for 2D rectangular visualization (XZ plane)
half_L = L_sat / 2;
half_H = H_sat / 2;
satellite_params.geometry.vertices_2d = [
    -half_L, -half_H;  % Bottom-left
     half_L, -half_H;  % Bottom-right
     half_L,  half_H;  % Top-right
    -half_L,  half_H   % Top-left
]';  % 2×4 matrix

%% Display Summary
fprintf('✓ Satellite 2D created from URDF:\n');
fprintf('  Mass: %.0f kg\n', robot_sat.base_link.mass);
fprintf('  Inertia Izz: %.1f kg·m²\n', robot_sat.base_link.inertia(3,3));
fprintf('  RW capacity: %.1f N·m·s (%.1f N·m max torque)\n', ...
        satellite_params.RW.h_max, satellite_params.RW.tau_max);
fprintf('  Contact points: %d docks at ±%.1f m\n', ...
        length(satellite_params.contact_points), 1.5);

end
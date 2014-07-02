%% HAWT Parameter file for CACTUS
% List of variables needed to generate a CACTUS *.geom file for a HAWT.

global CREATEGEOM_PATH

CREATEGEOM_PATH = '../../CreateGeom/';
addpath(CREATEGEOM_PATH)

%% Geometry
% Rotor Parameters
rotor_params.num_blades            = 3;				% Number of blades
rotor_params.radius                = 44.2913;			% Rotor radius (ft)
rotor_params.cone_angle            = 0;				% Cone angle (deg)
rotor_params.tilt_angle            = 0;				% Tilt angle (deg)

% Blade Parameters
blade_params.pitch                 = 0; % pitch (degrees) - positive is LE into the wind
blade_params.eta                   = 0; % Blade mount point ratio ((distance behind leading edge of the blade mount point) / (root chord)) 

% Discretization ("grid") Parameters
grid_params.node_distribution_type = 'sin';		% Distribution of grid points ([uniform],tanh,sin,custom)
grid_params.num_elems              = 20; % Number of line elements
grid_params.r_over_R_start         = 0.00; % r/R of innermost grid node
grid_params.r_over_R_end           = 1.00; % r/R of outermost grid node


%% Input
% Relative path to airfoil directory

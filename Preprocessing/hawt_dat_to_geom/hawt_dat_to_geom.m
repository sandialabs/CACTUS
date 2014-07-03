% HAWT_DAT_TO_GEOM : Reads in an aerodynamic schedule for a HAWT blade and creates a CACTUS-formatted .geom file.

function [Turbine] = hawt_dat_to_geom(dat_filename, geom_filename, rotor_params, blade_params, grid_params, varargin)

	%% Load aerodynamic schedule from DAT file, pack into a struct
	aero_schedule_nodes_raw = load(dat_filename);

	aero_schedule_nodes.r_over_R   = aero_schedule_nodes_raw(:,1);
	aero_schedule_nodes.c_over_R   = aero_schedule_nodes_raw(:,2);
	aero_schedule_nodes.beta       = aero_schedule_nodes_raw(:,3);
	aero_schedule_nodes.airfoil_id = aero_schedule_nodes_raw(:,4);

	%% Map the aerodynamic schedule to computational line elements with specified node locations
	r_over_R_range = grid_params.r_over_R_end - grid_params.r_over_R_start;

	switch grid_params.node_distribution_type
		case 'tanh'
			% Use a tanh distribution on both root and tip
			r_over_R_node_comp = (r_over_R_range)*(tanh(linspace(-pi/2,pi/2,grid_params.num_elems+1))/2 + 0.5) + grid_params.r_over_R_start; 
			
		case 'sin'
			% Use a sin distribution on both root and tip
			r_over_R_node_comp = (r_over_R_range)*(sin(linspace(-pi/2,pi/2,grid_params.num_elems+1))/2 + 0.5) + grid_params.r_over_R_start; 
			
		case 'uniform'
			% Use a uniform distribution across blade span
			r_over_R_node_comp = (r_over_R_range)*linspace(0,1,grid_params.num_elems+1) + grid_params.r_over_R_start;
			
		case 'custom'
			% Use a custom distribution specified in argument
			r_over_R_node_comp = varargin{1};
	end

	%% Call function to interpolate aero schedule data onto computational grid
	aero_schedule_elems  = distribute_nodes(aero_schedule_nodes, r_over_R_node_comp);

	%% Unpack and concatenate where necessary the 'mapped' aerodynamic schedule
	r_over_R_gridnodes = [aero_schedule_elems.r_over_R_a(1:end) aero_schedule_elems.r_over_R_b(end)];	% Blade element positions (nodes) (NBElem+1 elemnts)
	c_over_R_gridnodes = [aero_schedule_elems.c_over_R_a(1:end) aero_schedule_elems.c_over_R_b(end)];   % Blade chord / turbine radius (NBElem+1 elements ordered root to tip)
	beta_gridcenters   = [aero_schedule_elems.beta_a(1:end)         aero_schedule_elems.beta_b(end)];   % Blade station twist distribution in degrees (w.r.t. rotor disk plane, positive LE into wind (-x))
	af_id_gridcenters  = aero_schedule_elems.airfoil_id;												% Airfoil ID distribution
	
	%% Unpack parameters
	num_blades = rotor_params.num_blades;            % Number of blades
	R     	   = rotor_params.radius;                % Rotor radius (ft)
	cone_angle = rotor_params.cone_angle;            % Cone angle (deg)
	tilt_angle = rotor_params.tilt_angle;            % Tilt angle (deg)

	pitch      = blade_params.pitch;                 % pitch (degrees) - positive is LE into the wind
	eta        = blade_params.eta;                   % Blade mount point ratio ((distance behind leading edge of the blade mount point) / (root chord)) 
	
	turbine_type = 'HAWT'

	% Set some parameters
	reference_R_ratio = 1;                           % rotor radius / reference ratio -- set to 1 for simplicity
	hub_r             = grid_params.r_over_R_start;  % set hub radius as the radial position of the start of the innermost grid element

	%% Create turbine (using modified CreateGeom scripts)
	Turbine = CreateTurbine(num_blades,
	                        grid_params.num_elems,
	                        0,
	                        0,
	                        R,
	                        [],
	                        [],
	                        [],
	                        'HAWT',
	                        reference_R_ratio,
	                        hub_r,
	                        r_over_R_gridnodes,
	                        c_over_R_gridnodes,
	                        beta_gridcenters,
	                        af_id_gridcenters,
	                        pitch,
	                        eta,
	                        cone_angle,
	                        tilt_angle);

	%% Write geometry to file
	WriteTurbineGeom(geom_filename, Turbine)
end
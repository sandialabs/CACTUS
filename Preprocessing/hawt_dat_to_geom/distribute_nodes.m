% distribute_nodes.m
%
% DISTRIBUTE_NODES : 	Creates a 1-D distribution of nodes (r/R, c/R, beta, airfoil ID) along a blade surface.
% 						Linearly interpolates over an aerodynamic schedule which has data points defined at NODES.
%						Returns a data structure containing data for the discretized blade, indexed by blade element.
%
%						TL;DR - NODE Aerodynamic Schedule --> ELEMENT Aerodynamic Schedule with desired spacing

function [aero_schedule_elems] = distribute_nodes(aero_schedule_nodes, r_over_R_endpoints_elems)

	% Compute number of line elements
	num_elements = length(r_over_R_endpoints_elems)-1;

	% Get the locations of the element centers
	r_over_R_centers_elems = 0.5*(r_over_R_endpoints_elems(1:end-1) + r_over_R_endpoints_elems(2:end));

	
	% Break out data from structure
	r_over_R_nodes   = aero_schedule_nodes.r_over_R;
	c_over_R_nodes   = aero_schedule_nodes.c_over_R;
	beta_nodes       = aero_schedule_nodes.beta;
	airfoil_id_nodes = aero_schedule_nodes.airfoil_id;

	% Error Catching
	if length(r_over_R_endpoints_elems) < 2
		error('Too few nodes specified. At least one element (two or more endpoints) must be specified.')
	end


	% Warn about potential extrapolation
    if r_over_R_endpoints_elems(1) < r_over_R_nodes(1) || r_over_R_endpoints_elems(end) > r_over_R_nodes(end)
    	warning('Specified element endpoints are outside of given node data. Data will be linearly extrapolated.');
    end

	if r_over_R_centers_elems(1) < r_over_R_nodes(1) || r_over_R_centers_elems(end) > r_over_R_nodes(end)
		warning('Specified element endpoints produce element centers outside of given node data. Data will be linearly extrapolated.');
	end

	% Linearly interpolate to get data at desired points
	c_over_R_endpoints_elems = interp1(r_over_R_nodes, c_over_R_nodes, r_over_R_endpoints_elems, 'extrap');
	beta_endpoints_elems     = interp1(r_over_R_nodes,     beta_nodes, r_over_R_endpoints_elems, 'extrap');
	
	% Center values (unneeded, for now)
	c_over_R_centers_elems   = interp1(r_over_R_nodes, c_over_R_nodes, r_over_R_centers_elems, 'extrap');
	beta_centers_elems       = interp1(r_over_R_nodes,     beta_nodes, r_over_R_centers_elems, 'extrap');
	
	% For Airfoil ID, round to next lower node
	for i = 1:num_elements
		try 
			airfoil_id_elems(i) = airfoil_id_nodes(max(find(r_over_R_nodes < r_over_R_centers_elems(i))));
		catch
			warning(['No node location was found at a radius lower than the specified element r/R of ' num2str(r_over_R_centers_elems(i)) '. Defaulting to specified airfoil at innermost location.']);
			airfoil_id_elems(i) = airfoil_id_nodes(1);
		end
	end

	%% Pack into a structure (of arrays). Indexed by element number.
	aero_schedule_elems.num_elements    = length(r_over_R_centers_elems);
	aero_schedule_elems.r_over_R_center = r_over_R_centers_elems;
	aero_schedule_elems.r_over_R_a      = r_over_R_endpoints_elems(1:end-1);
	aero_schedule_elems.r_over_R_b      = r_over_R_endpoints_elems(2:end);
	aero_schedule_elems.beta_center     = beta_centers_elems;
	aero_schedule_elems.beta_a          = beta_endpoints_elems(1:end-1);
	aero_schedule_elems.beta_b          = beta_endpoints_elems(2:end);
	aero_schedule_elems.c_over_R_center = c_over_R_centers_elems;
	aero_schedule_elems.c_over_R_a      = c_over_R_endpoints_elems(1:end-1);
	aero_schedule_elems.c_over_R_b      = c_over_R_endpoints_elems(2:end);
	aero_schedule_elems.airfoil_id      = airfoil_id_elems;

end
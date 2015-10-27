% Author: Dexter Barrows
% McMaster University
% M745 - Topics in Numerical Analysis
% Final Project

% Implementation of Runge_Kutta 4 method for integration of a system of ODEs

% f 		- function handle for the rhs of the DE
% tspan 	- vector with start end end times for the integration
% y0 		- the initial condition
% h 		- step size

function y  = rk4(f,tspan,y0,h,pset,ah)

	t0 	= tspan(1);
	T 	= tspan(2);

	% make y0 a column vector if it is passed in as a row
	if size(y0, 1) ~= 1
   		y0 = y0'; 
	end

	% get size of system, number of steps required
	N = length(y0);
	nsteps = floor((T-t0)/h) + 1;

	% allocate storage for results, populate initial conditions

	t = t0;
	y = y0;
	i = 2;
	while t < T
		% update solution in time
		k1 = h * f(ah, y, pset);
		k2 = h * f(ah, y + k1/2, pset);
		k3 = h * f(ah, y + k2/2, pset);
		k4 = h * f(ah, y + k3, pset);
		y  = y + k1/6 + k2/3 + k3/3 + k4/6;

		t  = t + h;

		% increment
		i = i + 1;
	end

end
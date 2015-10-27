function PF

	nparticles = 300;
	data = textread('reg_1_data.txt');
	ah_dat = textread('nyc_qclim79_02.csv');

	% reconstruct weekly humidity forcing data
	ah_dat = ah_dat(round(linspace(1,365,52)));

	% particle states - each row is in format [S I], corresponds to a particular particle
	state(:,2) = 44*ones(nparticles,1); 	% infected
	state(:,1) = 14.6805e6 - state(:,2); 	% susceptible

	% base for perturbed parameters set in order [L D Rmax Rmin]
	pset = repmat([3.86 2.27 3.79 0.97], nparticles, 1);
	pert_vals = rand(nparticles,4) - 0.5;
	pset = pset + pert_vals;

	format long

	av_S = zeros(length(data));
	av_I = zeros(length(data));

	ah_counter = 40;

	len = 10

	all_states = zeros(nparticles,len);

	for i = 1:len

		obv = data(i);
		ah = ah_dat(ah_counter);
		ah_counter = mod(ah_counter+1,52) + 1;

		for k = 1:nparticles
			state(k,:) = rk4(@dSIR,[0 1/52],state(k,:),1/365,pset(k,:), ah);
		end

		av_S(i) = mean(state(:,1));
		av_I(i) = mean(state(:,2));
		all_states(:,i) = state(:,2);
		dev = std(state(:,2));

		weights = 1/(dev*sqrt(2*pi)) * exp( - (av_I(i) - state(:,2)).^2 / (2*dev^2) );
		weights = weights / sum(weights);

	end
	
	xx = 1:len;

	figure
	plot(xx,all_states);

	figure
	plot(xx(1:len),data(1:len),xx(1:len),av_I(1:len));
	legend('Data', 'PF');

	

end

function dX = dSIR(ah,y,pset)

	dX = zeros(1,2);
	S = y(1);
	I = y(2);

	L 		= pset(1);
	D 		= pset(2);
	Rmax    = pset(3);
	Rmin	= pset(4);

	a = -180;
	N = 14.6805e6;
	alpha = 0;
	b = log(Rmax - Rmin);
	Beta = 1*( (exp(a*ah + b) + Rmin) / D );

	dX(1) = (N-S-I)/L - Beta*I*S/N - alpha;
	dX(2) = Beta*I*S/N - I/D + alpha;

end
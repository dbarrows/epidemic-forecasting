function [xx,data,XX,proj_data] = Simplex2 (fac,shift)

	E = 5

	%data = textread('reg_1_data.txt');
	data = textread('Canada_gft_ili.dat');

	N = length(data);
	start_index = round(N*fac + shift);

	num_proj_steps = floor(start_index/2);

	proj_data = [data(1:start_index) ; zeros(num_proj_steps,1)];

	% construct library
	library = proj_data(E:start_index);
	for i = 1:(E-1)
		library = [ library proj_data(E-i:start_index-i) ];
	end

	% length of the library
	liblen = length(library);

	% determine largst value in time series for scaling against 52 weeks in a year
	max_ts_val = max(data);

	week_weights = mod( abs( (E:start_index) -  start_index ) , 52);
	week_weights(week_weights > 26) = 52 - week_weights(week_weights > 26);

	% scale week weightings so that they have as much impact as time series values
	week_weights = week_weights / 52 * max_ts_val;

	library = [ library week_weights' ];

	pre_x = proj_data(start_index);
	pre_vec = library(liblen,:);

	% compute norms and sort
	norms = sqrt( sum( (repmat(pre_vec,liblen,1) - library).^2 , 2) );
	[norm_asc, indices] = sort(norms);

	i = 2;
	nb_count = 0;
	nb_inds = zeros(E+1,1);
	while(nb_count < (E+1))
		if indices(i) < ( liblen - (num_proj_steps - 1) )
			nb_count = nb_count + 1;
			nb_inds(nb_count) = indices(i);
		end
		i = i + 1;
	end

	ref_norm = norms(nb_inds(1)); % get distance to nearest neighbor
	if ref_norm ~= 0
		weights = exp(-norms(nb_inds)/ref_norm);
	else
		weights = exp(-norms(nb_inds));
	end

	if (sum(weights) == 0)
		nb_inds
		ref_norm
		weights
		return
	end


	for p = 1:num_proj_steps

		x_next = sum(weights.*proj_data(nb_inds+E-1+p)) / sum(weights);

		proj_data(start_index+p) = x_next;

		if isnan(x_next)
			fprintf('Bad things')
			weights
			weights_sum = sum(weights)
			ref_norm
			nb_inds
			proj_data(nb_inds+E)
			return
		end

	end

	xx = 1:N;
	XX = start_index:(start_index+num_proj_steps);
	%figure
	%plot( xx,data,XX,proj_data(XX),'-', nb_inds+E,proj_data(nb_inds+E),'o' )
	proj_data= proj_data(XX);

	%hold all
	%scatter(library(:,2),library(:,1),10,'filled')
	%scatter(proj_data(nb_inds+E-1),proj_data(nb_inds+E),100,'filled')
	%scatter(pre_vec(2), pre_vec(1), 100,'filled')
	%scatter(nearby_vecs(:,2), nearby_vecs(:,1), 50, 'filled')
	%scatter(pre_x,x_next,200,'filled')


end
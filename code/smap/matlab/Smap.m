function Smap

	E = 4
	theta = 10

	data = textread('reg_1_data.txt');
	%data = textread('Canada_gft_ili.dat');

	N = length(data);
	start_index = N/2; %- 2*52;

	%num_proj_steps = min( N - start_index, floor(start_index/2) );
	num_proj_steps = 6*52

	proj_data = [data(1:start_index) ; zeros(num_proj_steps,1)];

	% construct library
	library = proj_data(E:start_index);
	for i = 1:(E-1)
		library = [ library proj_data(E-i:start_index-i) ];
	end

	% length of the library
	liblen = length(library);

	pre_x = proj_data(start_index);
	pre_vec = library(liblen,:);

	% compute norms
	norms = sqrt( sum( (repmat(pre_vec,liblen,1) - library).^2 , 2) );
	%[norm_asc, indices] = sort(norms)
	d_bar = mean(norms);

	for p = 1:num_proj_steps

		% max library index that can be included to make the prediction p steps ahead
		lib_max_index = liblen - p;

		lib_part = library(1:lib_max_index,:);
		norms_part = norms(1:lib_max_index);

		weights = exp(-theta*norms_part / d_bar);

		if (sum(weights) == 0)
			weights
			return
		end

		% indices that can be used in prediction
		nb_inds = 1:lib_max_index;

		% get b(i) = w(|x_i  - x_t|)*y_i
		b = weights .* proj_data(nb_inds+E-1+p);

		% get A(i,j) = w(|x_i  - x_t|)*x_i(j)
		A = repmat(weights,1,E) .* lib_part;

		c = A\b;

		x_next = sum(c .* pre_vec');

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
	figure
	plot( xx,data,XX,proj_data(XX),'-')

	%hold all
	%scatter(library(:,2),library(:,1),10,'filled')
	%scatter(proj_data(nb_inds+E-1),proj_data(nb_inds+E),100,'filled')
	%scatter(pre_vec(2), pre_vec(1), 100,'filled')
	%scatter(nearby_vecs(:,2), nearby_vecs(:,1), 50, 'filled')
	%scatter(pre_x,x_next,200,'filled')


end






















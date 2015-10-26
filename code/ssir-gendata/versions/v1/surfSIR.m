function surfSIR(dim)

	system('make');
	system( sprintf('./ssir %d', dim) );

	pause(2)

	data = textread('ssir_data.dat');

	[ydim,xdim] = size(data);

	xx = 1:xdim;

	[X,Y] = meshgrid(xx,xx);

	axis_bounds = [1 xdim 1 xdim 0 500];
	cmax = max(max(data));
	cmin = 0;

	figure

	Z = data(1:dim,:);

	surf(X,Y,Z);
		axis(axis_bounds);
		caxis([cmin cmax]);
	drawnow

	i = 1;
	while (1)

		ind = mod(i,(ydim/dim)-1) + 1;
		
		lbound = (ind-1)*dim+1;
		ubound = ind*dim;
		Z = data(lbound:ubound,:);

		surf(X,Y,Z);
			axis(axis_bounds);
			caxis([cmin cmax]);
		drawnow

		i = i + 1;

	end


end
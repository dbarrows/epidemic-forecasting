function surfSIR(filename)

	data = textread(filename);

	[ydim,xdim] = size(data);

	xx = 1:xdim;

	[X,Y] = meshgrid(xx,xx);

	axis_bounds = [1 xdim 1 xdim 0 500];
	cmax = max(max(data));
	cmin = 0;

	figure

	Z = data(1:xdim,:);

	surf(X,Y,Z);
		axis(axis_bounds);
		caxis([cmin cmax]);
		colormap(cool)
	drawnow

	i = 1;
	while (1)

		ind = mod(i,(ydim/xdim)-1) + 1;
		
		lbound = (ind-1)*xdim+1;
		ubound = ind*xdim;
		Z = data(lbound:ubound,:);

		surf(X,Y,Z);
			axis(axis_bounds);
			caxis([cmin cmax]);
			colormap(cool)
		drawnow

		i = i + 1;

	end


end
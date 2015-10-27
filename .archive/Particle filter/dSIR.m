function dX = dSIRflu(t,y,L,D,Rmax,Rmin)

	dX = zeros(1,2);
	S = y(1);
	I = y(2);

	N = 14.6805e6;
	alpha = 0;
	Beta = Rmax / D;

	dX(1) = (N-S-I)/L - Beta*I*S/N - alpha;
	dX(2) = Beta*I*S/N - I/D + alpha;

end
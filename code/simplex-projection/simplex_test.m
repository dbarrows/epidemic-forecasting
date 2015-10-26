fac = 1/2
shift = 52

[xx,data,XX,proj1] = Simplex(fac,shift);
[~,~,~,proj2] = Simplex2(fac,shift);

plot(xx,data,XX,proj1,XX,proj2)
legend('Data','Simplex','Simplex 2')
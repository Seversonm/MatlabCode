%
gamma = [-3:0.1:3];

figure
hold on

for dPrime = [0.5, 1.5, 3]
x = normcdf((-dPrime/2)- gamma);
y = normcdf((dPrime/2)- gamma);
plot(x,y)

end


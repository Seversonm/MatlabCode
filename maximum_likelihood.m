data(:,1) = [0,6.67,13.33,20,26.67,33.33,40,46.67]'; % x or abscissa values
data(:,2) = [0,3,10,34,57,67,71,72]'; % y or ordinate values
%starting parameter values
params0 = [20,4,1];

f = @(params)minimize(params);
params = fminsearch(f,params0);

%params = fmincon(@minimize,params0,[],[]);

prediction = psym(data,params);
plot(data(:,1), data(:,2)/72, 'o')
hold on
plot(data(:,1), prediction, 'r')

function [L] = minimize(params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% DON'T DO THIS EVER. IT BAD.
data(:,1) = [0,6.67,13.33,20,26.67,33.33,40,46.67]'; % x or abscissa values
data(:,2) = [0,3,10,34,57,67,71,72]'; % y or ordinate values
% BAD THINGS END HERE.
data(:,3) = 72 - data(:,2);

u = params(1);
a = params(2);
b = params(3);

L = -2*sum(data(:,2).*log(psym(data,params)) + data(:,3).*log(1-psym(data,params)));
end

function [F] = psym(data,params)
x = data(:,1);
u = params(1);
a = params(2);
b = params(3);

F = normcdf(0.5*sign(x - u).*((abs(x-u)/a).^b));
end

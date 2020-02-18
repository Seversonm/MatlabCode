figure(1)
x = (1:0.5:5);
y = (1:0.5:5);
plot(x,y,'b')

title ('Discriminability vs. target amplitude'); 
xlabel ('target amplitude');
ylabel ('dPrime');

figure(2)
x = (1:0.5:5);
y = (1:0.5:5);
plot(x,y,'r')

title ('Amplitude threshold vs. noise standard deviation'); 
xlabel ('sigma');
ylabel ('alpha threshold');


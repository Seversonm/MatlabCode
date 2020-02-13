p.noiseMean = 0;   
p.signalMean = 1;
p.sd = 1;
z = -4:.2:6;  %response-axis sampling points

noise.y  = normpdf(z,p.noiseMean,p.sd);
signal.y = normpdf(z,p.signalMean,p.sd);

ylim = get(gca,'YLim');

p.criterion = 1.5;

pHit = 1-normcdf(p.criterion,p.signalMean,p.sd)
pFA =  1-normcdf(p.criterion,p.noiseMean,p.sd)

dPrime = 0.5       %(p.signalMean-p.noiseMean)/p.sd; Given in this problem. 
fprintf(sprintf('  dPrime =    %5.2f\n',dPrime));

pHits = 1-normcdf(z,p.signalMean,p.sd);
pFAs  = 1-normcdf(z,p.noiseMean,p.sd);
figure(1)
hold on
axis equal
axis tight
xlabel('pFA');
ylabel('pHit');
plot(pFAs,pHits,'k-');

hold on
dPrime = 1.5;
plot(pFAs,pHits,'r-');

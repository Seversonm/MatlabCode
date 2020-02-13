p.noiseMean = 0;   
p.signalMean = 1;
p.sd = 1;
z = -4:.2:6;
p.criterion = 1.5;

pHits = 1-normcdf(z,p.signalMean,p.sd);
pFAs  = 1-normcdf(z,p.noiseMean,p.sd);
figure(1)
clf
hold on
axis equal
axis tight
xlabel('pFA');
ylabel('pHit');
plot(pFAs,pHits, 'k-');
hold on
p.criterion = 2;
plot(pFAs,pHits, 'r-');
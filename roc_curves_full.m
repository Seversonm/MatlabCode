p.noiseMean = 0;   
p.signalMean = 1;
p.sd = 1;
z = -4:.2:6;  %response-axis sampling points

noise.y  = normpdf(z,p.noiseMean,p.sd);
signal.y = normpdf(z,p.signalMean,p.sd);

figure(1)
clf
plot(z,noise.y);
hold on
plot(z,signal.y,'r-');

ylim = get(gca,'YLim');

text(p.noiseMean,ylim(2)*.9,'Noise','VerticalAlignment','top','HorizontalAlignment','center','Color','b');
text(p.signalMean,ylim(2)*.9,'Signal','VerticalAlignment','top','HorizontalAlignment','center','Color','r');
xlabel('Internal Response');
p.criterion = 1.5;
plot(p.criterion*[1,1],ylim,'k:');

hcrit = text(p.criterion,0,'criterion','VerticalAlignment','bottom','HorizontalAlignment','left');
pHit = 1-normcdf(p.criterion,p.signalMean,p.sd)
pFA =  1-normcdf(p.criterion,p.noiseMean,p.sd)
disp(' ');
disp('           |         Response     |')
disp('  Signal   |  "Yes"    |  "No"    |')
disp('  ---------------------------------')
fprintf('  Present  |   %3.1f%%   |   %3.1f   |\n',100*pHit,100*(1-pHit));
disp('  ---------+-----------+----------|');
fprintf('  Absent   |   %3.1f    |   %3.1f   |\n',100*pFA,100*(1-pFA));
disp('  ---------------------------------');
PC = (pHit + (1-pFA))/2;  %proportion correct
fprintf('  Percent Correct: %5.2f%%\n',100*PC);
dPrime = 0.5       %(p.signalMean-p.noiseMean)/p.sd; Given in this problem. 
fprintf(sprintf('  dPrime =    %5.2f\n',dPrime));

pHits = 1-normcdf(z,p.signalMean,p.sd);
pFAs  = 1-normcdf(z,p.noiseMean,p.sd);
figure(2)
clf
hold on
axis equal
axis tight
xlabel('pFA');
ylabel('pHit');
plot(pFAs,pHits,'k-');
% Retinotopy - works best when you image the vertical bars first and then
% the horizontal bars. This version is for when you split the recording between the vertical bars 
% and the horizontal bars, and this particular analysis is for the VERTICAL bars. 
% First open your data in  Caliman and export, then run the following code.
% Tested 1/30/2020. 

for n=1:size(data.vdat,2);vsig(n)=mean(data.vdat{n});end;

detect = 0;
x=1;
for n=1:size(vsig,2);
   
    if vsig(n) > 2000 & detect == 0;
        if size(find(vsig(1:n)<1500),1)>0;
        vsig_times(x) = max(find(vsig(1:n)<1500));
        else
        vsig_times(x) = 0;
        end;
        detect = 1;
        x = x+1;
    end;
    if vsig(n) < 1500 & detect == 1;
        detect = 0;
    end;

end;
if detect == 0 & x == 1;
    vsigtimes=[];
end;

%%

FR = 7910/256;
for n=1:10
    k = (n-1)*6+1;
    onset_times(k)=vsig_times(n)+FR/2;
    for m=1:5
        onset_times(k+m)=onset_times(k+m-1)+6*FR;
    end;
end;
stim_time = round(onset_times);

I(1:256,1:455,1:size(vsig,2))=0;
for n=1:size(vsig,2);I(:,:,n)=data.green{n};end;
Ineuro(1:256,1:455,1:124,1:115)=0;
for m=1:size(stim_time,2);
    disp(m);
    Ineuro(:,:,1:124,m) = I(:,:,stim_time(m)-62:stim_time(m)+61);
end;

spatfltr = ones(50,50)/2500;
imagesc(imfilter(mean(mean(Ineuro(:,:,63:124,:),4),3),spatfltr)-imfilter(mean(mean(Ineuro(:,:,1:62,:),4),3),spatfltr))

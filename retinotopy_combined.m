% Retinotopy - works best when you image the vertical bars first and then
% the horizontal bars. This version is for when you record the vertical bars 
% and then horizontal bars in the same file.
% First open your data in  Caliman and export, then run the following code.

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
T=k+m+1;

for n=11:20
    k = T+(n-11)*5;
    onset_times(k)=vsig_times(n)+FR/2;
    for m=1:4
        onset_times(k+m)=onset_times(k+m-1)+6*FR;
    end;
end;
T=k+m+1;
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

% Modify the last line to display individual conditions (current format is to produce an activity map averaged across all conditions).


%% To see what the activity map looks like at every position where the bars were, use the following code: 
% Note: This example only shows six positions because it was run on a file in which
% only the vertical bars were recorded. If both vertical and horizontal were recorded in the same file, 
% you should have 11 positions (6 vertical and 5 horizontal).

% imagesc(imfilter(mean(mean(Ineuro(:,:,63:124,1:60),4),3),spatfltr)-imfilter(mean(mean(Ineuro(:,:,1:62,1:60),4),3),spatfltr))

% This is for position 1 of 6. If there are 11 positions, just change the line below to say "1:11:110"
% imagesc(imfilter(mean(mean(Ineuro(:,:,63:124,1:6:60),4),3),spatfltr)-imfilter(mean(mean(Ineuro(:,:,1:62,1:6:60),4),3),spatfltr))
% figure 
% imagesc(imfilter(mean(mean(Ineuro(:,:,63:124,2:6:60),4),3),spatfltr)-imfilter(mean(mean(Ineuro(:,:,1:62,2:6:60),4),3),spatfltr))
% figure
% imagesc(imfilter(mean(mean(Ineuro(:,:,63:124,3:6:60),4),3),spatfltr)-imfilter(mean(mean(Ineuro(:,:,1:62,3:6:60),4),3),spatfltr))
% figure
% imagesc(imfilter(mean(mean(Ineuro(:,:,63:124,4:6:60),4),3),spatfltr)-imfilter(mean(mean(Ineuro(:,:,1:62,4:6:60),4),3),spatfltr))
% figure
% imagesc(imfilter(mean(mean(Ineuro(:,:,63:124,5:6:60),4),3),spatfltr)-imfilter(mean(mean(Ineuro(:,:,1:62,5:6:60),4),3),spatfltr))
% figure
% imagesc(imfilter(mean(mean(Ineuro(:,:,63:124,6:6:60),4),3),spatfltr)-imfilter(mean(mean(Ineuro(:,:,1:62,6:6:60),4),3),spatfltr))
% figure

% The line below takes out the filter so you can see the vasculature.
% imagesc(mean(mean(Ineuro(:,:,1:124,1:60),4),3)) 

% The line below does a left/right flip to get it matched with looking down
% on the window. The bottom of the image is lateral and the top is medial. 
% imagesc(fliplr(imfilter(mean(mean(Ineuro(:,:,63:124,4:6:60),4),3),spatfltr)-imfilter(mean(mean(Ineuro(:,:,1:62,4:6:60),4),3),spatfltr)))
% imagesc(fliplr(mean(mean(Ineuro(:,:,1:124,1:60),4),3)))

























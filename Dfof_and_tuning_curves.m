st = size(unique(data.stim(:,2)'),2); %% finds the number of unique stimuli played in this program. 
for i = 1:st %% for each unique stimuli
    inds(:,i) = find(data.stim(:,2)==i);%% creates an array. Finds the positions in the original data.stim file where this stimulus was played.
    frames(:,i) = data.stim(inds(:,i),1);%% creates an array. Finds the values (frame in which stimulus started) that go with each position in inds.
end
nt = length(frames(:,1)); %% gives the length of the frames array (number of repeats).
cellnum = length(data.Fraw(:,1)); %% Tells us how many rows of data there are in data.Fraw. This corresponds to how many cells were circled. 
for s = 1:st %% make a matrix with stim number, repeats, and time. This means for each unique stimuli.
    for n = 1:nt %% and for each frame in which the stimulus started
        for c = 1:cellnum %% and for each cell marked as a region of interest
            response(c,s,n,:) = squeeze(data.Fraw(c,frames(n,s):frames(n,s)+192)); %% find the response and create an array with stim number, repeats, and time. 
        end
    end
end

%%%%%%STUFF YOU NEED TO EDIT%%%%%%%%
cellnum1 = cellnum-1; %we don't care about the blank
elements = data.stim(:,2);
elements = elements'; 
%%%%%%%%%%Set some times%%%%%%%%%%%%
delay = 3; % time from ttl->show. What is the delay before the stim shows? 3 sec for my protocol.
bas = delay+32; %for delta F
stimend = 352; %end of stim + some more
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%preallocate variables%%%%%%%
combFraw = zeros(nt,((stimend-delay)+1),st);
combf = zeros(nt,((stimend-delay)+1),st);
degmean = zeros(st,1);
degstd = zeros(st,1);
degstderr = zeros(st,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for b = 1:cellnum1 %cell number. set the initial conditions equal to 1 to test. 
    w = cellnum; % blank dark cell
    n = nt;
for p = 1:n %repeat number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%raw signal for each trial%%%%%    
for i = 1:st
    combFraw(p,:,i) = data.Fraw(b,frames(p,i)+delay:frames(p,i)+stimend)-data.Fraw(w,frames(p,i)+delay:frames(p,i)+stimend);                   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%delta F calc%%%%%%%%%%%%%%
for p=1:n
     for i=1:st
combf(p,:,i)=(combFraw(p,:,i)-mean(combFraw(p,delay:bas,i)))./mean(combFraw(p,delay:bas,i));
     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%frames to seconds%%%%%%%%%%%%
time(1) = 1/32; %32 frames/s
for i =2:((stimend-delay)+1) % 2 to end of frame
    time(i) = time(i-1) + 1/32;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%create trial averages%%%%%%%%%%
combfa = mean(combf,1);
combFrawa = mean(combFraw,1);
combstd = std(combFraw,1);
combstderr = combstd/(sqrt(cellnum1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%Plot Raw data%%%%%%%%%%%%%% 
 figure (b+200)
 for j=1:st %num of stim. condition
 subplot(st,1,j)
    for i=1:n          
 plot(time,combFraw(i,:,j),'k')
    hold on
 plot(time,combFrawa(:,:,j),'r')
%     ylim([-1 1])
     xlim([0 11])
     line([1.5 1.5],[0 200]) %stim onset
    hold on
     line([5.5 5.5],[0 200])
set(gcf, 'Position', [100, 100, 400, 850])
    end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%Plot Delta F%%%%%%%%%%%%%%
figure (b+300)
  for j=1:st %num of stim. condition
  H(j)=subplot(st,1,j); %H contains axis handle
  hold on
     for i=1:n          
  plot(time,combf(i,:,j),'k')
  hold on
  plot(time,combfa(:,:,j),'r')
  hold on
 %ylim([-20 20])
xlim([0 11])
%line([1.5 1.5],[-1 1]) %stim onset
 % hold on
%line([5.5 5.5],[-1 1]) %stim offset
set(gcf, 'Position', [100, 100, 400, 850])
 end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%plot tuning curve%%%%%%%%%%
%%Tuning curve from raw signal
for i = 1:st
    degmean(i,1) = mean(combFrawa(1,128:256,i));
end
degmean = degmean';
degmeanf = degmean(1,1:8);
for i = 1:st
     degstderr(i,1) = mean(combstd(1,128:256,i));
end 
degstderr = degstderr';
degstderrf = degstderr(1,1:8);
ax = [0,45,90,135,180,225,270,315];
xlabels = {'0°','45°','90°','135°','180°','225°','270°','315°'};
figure(b+400)   
     errorbar(ax,degmeanf,degstderrf)
     set(gca, 'Xtick', ax, 'XtickLabel', xlabels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%Clean Workspace%%%%%%%%%%%%
%clear a ax b bas bl ce degmean degstd del el en i l n p w st rp xlabels y z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% F(t) = F(t) - Fo /Fo
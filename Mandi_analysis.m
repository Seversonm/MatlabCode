st = size(unique(data.stim(:,2)'),2); %% finds the number of unique stimuli played in this program. 
for i = 1:st %% for each unique stimuli
    inds(:,i) = find(data.stim(:,2)==i);%% creates an array. Finds the positions in the original data.stim file where this stimulus was played.
    frames(:,i) = data.stim(inds(:,i),1);%% creates an array. Finds the values (frame in which stimulus started) that go with each position in inds.
end
nt = length(frames(:,1)); %% gives the length of the frames array.
cellnum = length(data.Fraw(:,1)); %% Tells us how many rows of data there are in data.Fraw. This corresponds to how many cells were circled. 
for s = 1:st %% make a matrix with stim number, repeats, and time. This means for each unique stimuli.
    for n = 1:nt %% and for each frame in which the stimulus started
        for c = 1:cellnum %% and for each cell marked as a region of interest
            response(c,s,n,:) = squeeze(data.Fraw(c,frames(n,s):frames(n,s)+192)); %% find the response and create an array with stim number, repeats, and time. 
        end
    end
end
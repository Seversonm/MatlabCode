function [stim,gcstim, frcode] = getstim(gccode)
  debug = 1;
  stim = []; gcstim = []; frcode = [];
  
  str={'None','Flash Adaptation','Disparity','Drift Ori','SF','Flash Ori','oriattwosf', 'oldori', 'MT Drfting direction', 'MT Plaid 120', 'MT Plaid 90', 'MT Plaid 30 type II', 'MT SF tuning',  'MT TF tuning', 'Marmie'};
  [s,v] = listdlg('PromptString','StimulusType','SelectionMode','single','ListString',str);

  if (v==1) & (s>1)
        if s == 2
            [stim] = getadapstim();
        elseif s==3
            [stim] = getBINOCnewstim();
        elseif s==4
            [stim] = getDRORInewstim();
        elseif s==5
            [stim] = getSFnewstim();
        elseif s==6
            [stim] = getDRORInewstim();
        elseif s==7
            [stim] = getoriattwosf();
        elseif s==8
            [stim] = getORIstim2();
        elseif s==9
            [stim] = getMTdrdirectionstim();
        elseif s==10
            [stim] = getMTplaid120stim();
        elseif s== 11
            [stim] = getMTplaid90stim();
        elseif s==12
            [stim] = getMTplaidtypeIIstim();
        elseif s==13
            [stim] = getMTSFtuning();
        elseif s==14
            [stim] = getMTTFtuning();
        elseif s==15
            [stim] = getMarmieStim();

        end
      
      
        if ~isempty(stim)
            if debug
                figure(99)
                hist(stim(:,2),0:10); title(['stim numbers sorted  number of trials:',num2str(length(stim)/9)])
%                 figure(101)
%                 plot(x)
%                 hold on
%                 plot(stim(:,1),x(stim(:,1)),'*r')
%                 hold off
            end

            if ~isempty(gccode)
              % Find the frames that have the gccode
              frcode = [];
              for j = 1:length(gccode)
                if ~isempty(gccode{j})
                    frcode = [frcode j];
                end
              end
              for j = 1:length(frcode)
                gcstim{j} = gccode{frcode(j)};
              end
    
            end
        else
            stim=zeros(2,2);
        end

  end
return



      
function [stim] = getadapstim(gccode)
debug = 1;
% condition:
%   0   400000
%   1   402000
%   2   402200
%   3   402220
%   4   402020
%   5   400200
%   6   400220
%   7   400020
%   8   404000
%   9   404400
%   10  404440


global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

z = x;
x = diff(x);

x = resamplevec(x,samp);
z = resamplevec(z,samp);
x(x<0) = 0;
z(z<0) = 0;
x = x./max(x);
z = z./max(z);
z = round(z * 4);

nx = length(x);
rawpeak = (find((x(1:nx-2)<x(2:nx-1))&(x(2:nx-1)>x(3:nx)))+1)';
%drop things close to 0
rawpeak = rawpeak(x(rawpeak)>.7);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 

    
    % Let's read out the code for A, B, C, D
    a = 4;
    b= mean(z(onsets(ii)+(17:18)));
    c = mean(z(onsets(ii)+(23:24)));
    d = mean(z(onsets(ii)+(32:33)));

%    disp([onsets(ii) a b c d])

    e = 100*b + 10*c + d;
    if e==0
        stim = [stim;onsets(ii) 0];
    elseif e == 200
        stim = [stim; onsets(ii) 1];
    elseif e == 220
        stim = [stim; onsets(ii) 2];
    elseif e == 222 
        stim = [stim; onsets(ii) 3];
    elseif e == 202
        stim = [stim; onsets(ii) 4];
    elseif e == 20
        stim = [stim; onsets(ii) 5];
    elseif e == 22
        stim = [stim; onsets(ii) 6];
    elseif e == 2
        stim = [stim; onsets(ii) 7];
    elseif e == 400
        stim = [stim; onsets(ii) 8];
    elseif e == 440
        stim = [stim; onsets(ii) 9];
    elseif e == 444
        stim = [stim; onsets(ii) 10];
    else
        disp('I have no idea what that code is for!!!!')
    end
    
end



function [stim] = getBINOCnewstim()

global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

z = x;
x = diff(x);

x = resamplevec(x,samp);
z = resamplevec(z,samp);
x(x<0) = 0;
z(z<0) = 0;
x = x./max(x);
z = z./max(z);
z = round(z * 5);
x = round(x * 5);
nx = length(x);

rawpeak = (find((x(1:nx-2)<x(2:nx-1))&(x(2:nx-1)>x(3:nx)))+1)';
%drop things close to 0
%rawpeak = rawpeak(x(rawpeak)>.7);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
v = [];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 
    
    
    % Let's read out the code for A, B, C, D
    a = 5;
    b= mean(z(onsets(ii)+(12:13)));
    c = mean(z(onsets(ii)+(20:21)));
    d = mean(z(onsets(ii)+(29:30)));
    e = mean(z(onsets(ii)+(38:39)));

%    disp([onsets(ii) a b c d])

    val = 10000*a + 1000*b + 100*c + 10*d + e;
    v = [v; val a b c d e];
    if val==55000
        stim = [stim;onsets(ii) 1];
    elseif val == 55500
        stim = [stim; onsets(ii) 2];
    elseif val == 55550
        stim = [stim; onsets(ii) 3];
    elseif val == 55555 
        stim = [stim; onsets(ii) 4];
    elseif val == 50555
        stim = [stim; onsets(ii) 5];
    elseif val == 55055
        stim = [stim; onsets(ii) 6];
    elseif val == 55505
        stim = [stim; onsets(ii) 7];
    elseif val == 50505
        stim = [stim; onsets(ii) 8];
    elseif val == 50000
        stim = [stim; onsets(ii) 9];
    elseif val == 50005
        stim = [stim; onsets(ii) 10];
    elseif val == 50050
        stim = [stim; onsets(ii) 11];
    else
        disp('I have no idea what that code is for!!!!')
    end
    
end



function [stim] = getDRORInewstim()

global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

z = x;
x = diff(x);

x = resamplevec(x,samp);
z = resamplevec(z,samp);
x(x<0) = 0;
z(z<0) = 0;
x = x./max(x);
z = z./max(z);
z = round(z * 5);
x = round(x * 5);
nx = length(x);

rawpeak = (find((x(1:nx-2)<x(2:nx-1))&(x(2:nx-1)>x(3:nx)))+1)';
%drop things close to 0
%rawpeak = rawpeak(x(rawpeak)>.7);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
v = [];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 
    
    
    % Let's read out the code for A, B, C, D
    a = 5;
    b= mean(z(onsets(ii)+(12:13)));
    c = mean(z(onsets(ii)+(20:21)));
    d = mean(z(onsets(ii)+(29:30)));
    e = mean(z(onsets(ii)+(38:39)));

%    disp([onsets(ii) a b c d])

    val = 10000*a + 1000*b + 100*c + 10*d + e;
    v = [v; val a b c d e];
    if val==55000
        stim = [stim;onsets(ii) 1];
    elseif val == 55500
        stim = [stim; onsets(ii) 2];
    elseif val == 55550
        stim = [stim; onsets(ii) 3];
    elseif val == 55555 
        stim = [stim; onsets(ii) 4];
    elseif val == 50555
        stim = [stim; onsets(ii) 5];
    elseif val == 55055
        stim = [stim; onsets(ii) 6];
    elseif val == 55505
        stim = [stim; onsets(ii) 7];
    elseif val == 50505
        stim = [stim; onsets(ii) 8];
    elseif val == 50000
        stim = [stim; onsets(ii) 9];
    else
        disp('I have no idea what that code is for!!!!')
    end
    
end

function [stim] = getSFnewstim()

global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

z = x;
x = diff(x);

x = resamplevec(x,samp);
z = resamplevec(z,samp);
x(x<0) = 0;
z(z<0) = 0;
x = x./max(x);
z = z./max(z);
z = round(z * 5);
x = round(x * 5);
nx = length(x);

rawpeak = (find((x(1:nx-2)<x(2:nx-1))&(x(2:nx-1)>x(3:nx)))+1)';
%drop things close to 0
%rawpeak = rawpeak(x(rawpeak)>.7);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
v = [];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 
    
    
    % Let's read out the code for A, B, C, D
    a = 5;
    b= mean(z(onsets(ii)+(12:13)));
    c = mean(z(onsets(ii)+(20:21)));
    d = mean(z(onsets(ii)+(29:30)));
    e = mean(z(onsets(ii)+(38:39)));

%    disp([onsets(ii) a b c d])

    val = 10000*a + 1000*b + 100*c + 10*d + e;
    v = [v; val a b c d e];
    if val==55000
        stim = [stim;onsets(ii) 1];
    elseif val == 55500
        stim = [stim; onsets(ii) 2];
    elseif val == 55550
        stim = [stim; onsets(ii) 3];
    elseif val == 55555 
        stim = [stim; onsets(ii) 4];
    elseif val == 50555
        stim = [stim; onsets(ii) 5];
    elseif val == 55055
        stim = [stim; onsets(ii) 6];
    elseif val == 55505
        stim = [stim; onsets(ii) 7];
    elseif val == 50000
        stim = [stim; onsets(ii) 8];
    else
        disp('I have no idea what that code is for!!!!')
    end
    
end

function [stim] = getoriattwosf()

global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

z = x;
x = diff(x);

x = resamplevec(x,samp);
z = resamplevec(z,samp);
x(x<0) = 0;
z(z<0) = 0;
x = x./max(x);
z = z./max(z);
z = round(z * 5);
x = round(x * 5);
nx = length(x);
rawpeak = find((z(2:nx-1)==1)&(z(1:nx-2) == 0) & z(3:nx) == 1);
%drop things close to 0
%rawpeak = rawpeak(x(rawpeak)>.7);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
v = [];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 
    
    
    % Let's read out the code for A, B, C, D
    a = 1;
    b= mean(z(onsets(ii)+(8)));
    c = mean(z(onsets(ii)+(14:15)));
    d = mean(z(onsets(ii)+(20:21)));
    e = mean(z(onsets(ii)+(40:41)));

%    disp([onsets(ii) a b c d])

    val = 10000*a + 1000*b + 100*c + 10*d + e;
    v = [v; val a b c d e];
    if val==14440
        stim = [stim;onsets(ii) 2];
    elseif val == 12240
        stim = [stim; onsets(ii) 3];
    elseif val == 14400
        stim = [stim; onsets(ii) 4];
    elseif val == 12200 
        stim = [stim; onsets(ii) 5];
    elseif val == 14240
        stim = [stim; onsets(ii) 6];
    elseif val == 12040
        stim = [stim; onsets(ii) 7];
    elseif val == 14200
        stim = [stim; onsets(ii) 8];
    elseif val == 12000
        stim = [stim; onsets(ii) 9];
    elseif val == 14040
        stim = [stim; onsets(ii) 10];
    elseif val == 10440
        stim = [stim; onsets(ii) 11];
    elseif val == 14000
        stim = [stim; onsets(ii) 12];
    elseif val == 10400
        stim = [stim; onsets(ii) 13];
    elseif val == 10000
        stim = [stim; onsets(ii) 1];
    else
        disp('I have no idea what that code is for!!!!')
    end
    
end

if ~isempty(stim)
figure(99)
hist(stim(:,2),-1000:150); title(['stim numbers sorted  number of trials:',num2str(length(stim)/9)])
figure(101)
plot(x)
hold on
plot(stim(:,1),x(stim(:,1)),'*r')
else
    stim=zeros(2,2);
end

function [stim] = getORIstim2()

% currently good for boht orientaiton (9 stimuli) and counterphase gratings
% which is same # (1 blank, 4 phases at 0, 4 phases at 90)
%

% key:
%   1       blank
%   11      0 deg
%   111     45 deg
%   1111    90 deg
%   11111   135 deg
%   10111   180 deg
%   11011   225 deg
%   11101   270 deg
%   10101   315 deg
%   stim onset-offset ~40-50 frames
%   each pulse ~6-8 isi

global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

x = diff(x);
x = resamplevec(x,samp);
x(x<0) = 0;
x = x./max(x);
nx = length(x);
rawpeak = (find((x(1:nx-2)<x(2:nx-1))&(x(2:nx-1)>x(3:nx)))+1)';
%drop things close to 0
rawpeak = rawpeak(x(rawpeak)>.1);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 
    
    %TTL pulses go duration of imaging?
    if onsets(ii)+40 <= length(x)
        
    pulse = x(onsets(ii)-1:onsets(ii)+40);
    pulse = find(pulse>.95);
    
    if length(pulse)==1
        stim=[stim; onsets(ii)  1];
    elseif length(pulse)==2 
        stim=[stim; onsets(ii)  2];
    elseif length(pulse)==3 && sum(pulse(2:end)-pulse(1:end-1)<ipi)==2
        stim=[stim; onsets(ii)  3];
    elseif length(pulse)==4 && sum(pulse(2:end)-pulse(1:end-1)<ipi)==3
        stim=[stim; onsets(ii)  4];
    elseif length(pulse)==5
        stim=[stim; onsets(ii)  5];
    elseif length(pulse)==4 && pulse(2)-pulse(1)>=ipi
        stim=[stim; onsets(ii)  6];
    elseif length(pulse)==4 && pulse(3)-pulse(2)>=ipi
        stim=[stim; onsets(ii)  7];
    elseif length(pulse)==4 && pulse(4)-pulse(3)>=ipi
        stim=[stim; onsets(ii)  8]; 
        
    elseif length(pulse)==3 && sum(pulse(2:end)-pulse(1:end-1)>=ipi)==2
        stim=[stim; onsets(ii)  9];      
    end
    
    else
        disp 'dropped pulse'
    end
end

if ~isempty(stim)
figure(99)
hist(stim(:,2),1:9); title(['stim numbers sorted  number of trials:',num2str(length(stim)/9)])
figure(101)
plot(x)
hold on
plot(stim(:,1),x(stim(:,1)),'*r')
else
    stim=zeros(2,2);
end

function [stim] = getMTdrdirectionstim();

global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

z = x;
x = diff(x);

x = resamplevec(x,samp);
z = resamplevec(z,samp);
x(x<0) = 0;
z(z<0) = 0;
x = x./max(x);
z = z./max(z);
z = round(z * 5);
x = round(x * 5);
nx = length(x);
rawpeak = find((z(2:nx-1)==1)&(z(1:nx-2) == 0) & z(3:nx) == 1);
%drop things close to 0
%rawpeak = rawpeak(x(rawpeak)>.7);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
v = [];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 
    
    
    % Let's read out the code for A, B, C, D
    a = 1;
    b= mean(z(onsets(ii)+(8)));
    c = mean(z(onsets(ii)+(14:15)));
    d = mean(z(onsets(ii)+(20:21)));
    e = mean(z(onsets(ii)+(40:41)));

%    disp([onsets(ii) a b c d])

    val = 10000*a + 1000*b + 100*c + 10*d + e;
    v = [v; val a b c d e];
    if val==14440
        stim = [stim;onsets(ii) 2];
    elseif val == 12240
        stim = [stim; onsets(ii) 3];
    elseif val == 14400
        stim = [stim; onsets(ii) 4];
    elseif val == 12200 
        stim = [stim; onsets(ii) 5];
    elseif val == 14240
        stim = [stim; onsets(ii) 6];
    elseif val == 12040
        stim = [stim; onsets(ii) 7];
    elseif val == 14200
        stim = [stim; onsets(ii) 8];
    elseif val == 12000
        stim = [stim; onsets(ii) 9];
    elseif val == 14040
        stim = [stim; onsets(ii) 10];
    elseif val == 10440
        stim = [stim; onsets(ii) 11];
    elseif val == 14000
        stim = [stim; onsets(ii) 12];
    elseif val == 10400
        stim = [stim; onsets(ii) 13];
    elseif val == 10000
        stim = [stim; onsets(ii) 1];
    else
        disp('I have no idea what that code is for!!!!')
    end
    
end

if ~isempty(stim)
figure(99)
hist(stim(:,2),-1000:150); title(['stim numbers sorted  number of trials:',num2str(length(stim)/9)])
figure(101)
plot(x)
hold on
plot(stim(:,1),x(stim(:,1)),'*r')
else
    stim=zeros(2,2);
end

function [stim] = getMTplaid120stim();

global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

z = x;
x = diff(x);

x = resamplevec(x,samp);
z = resamplevec(z,samp);
x(x<0) = 0;
z(z<0) = 0;
x = x./max(x);
z = z./max(z);
z = round(z * 5);
x = round(x * 5);
nx = length(x);
rawpeak = find((z(2:nx-1)==1)&(z(1:nx-2) == 0) & z(3:nx) == 1);
%drop things close to 0
%rawpeak = rawpeak(x(rawpeak)>.7);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
v = [];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 
    
    
    % Let's read out the code for A, B, C, D
    a = 1;
    b= mean(z(onsets(ii)+(8)));
    c = mean(z(onsets(ii)+(14:15)));
    d = mean(z(onsets(ii)+(20:21)));
    e = mean(z(onsets(ii)+(40:41)));

%    disp([onsets(ii) a b c d])

    val = 10000*a + 1000*b + 100*c + 10*d + e;
    v = [v; val a b c d e];
    if val==14440
        stim = [stim;onsets(ii) 2];
    elseif val == 12240
        stim = [stim; onsets(ii) 3];
    elseif val == 14400
        stim = [stim; onsets(ii) 4];
    elseif val == 12200 
        stim = [stim; onsets(ii) 5];
    elseif val == 14240
        stim = [stim; onsets(ii) 6];
    elseif val == 12040
        stim = [stim; onsets(ii) 7];
    elseif val == 14200
        stim = [stim; onsets(ii) 8];
    elseif val == 12000
        stim = [stim; onsets(ii) 9];
    elseif val == 14040
        stim = [stim; onsets(ii) 10];
    elseif val == 10440
        stim = [stim; onsets(ii) 11];
    elseif val == 14000
        stim = [stim; onsets(ii) 12];
    elseif val == 10400
        stim = [stim; onsets(ii) 13];
    elseif val == 10000
        stim = [stim; onsets(ii) 1];
    else
        disp('I have no idea what that code is for!!!!')
    end
    
end

if ~isempty(stim)
figure(99)
hist(stim(:,2),-1000:150); title(['stim numbers sorted  number of trials:',num2str(length(stim)/9)])
figure(101)
plot(x)
hold on
plot(stim(:,1),x(stim(:,1)),'*r')
else
    stim=zeros(2,2);
end

function [stim] = getMTplaid90stim();

global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

z = x;
x = diff(x);

x = resamplevec(x,samp);
z = resamplevec(z,samp);
x(x<0) = 0;
z(z<0) = 0;
x = x./max(x);
z = z./max(z);
z = round(z * 5);
x = round(x * 5);
nx = length(x);
rawpeak = find((z(2:nx-1)==1)&(z(1:nx-2) == 0) & z(3:nx) == 1);
%drop things close to 0
%rawpeak = rawpeak(x(rawpeak)>.7);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
v = [];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 
    
    
    % Let's read out the code for A, B, C, D
    a = 1;
    b= mean(z(onsets(ii)+(8)));
    c = mean(z(onsets(ii)+(14:15)));
    d = mean(z(onsets(ii)+(20:21)));
    e = mean(z(onsets(ii)+(40:41)));

%    disp([onsets(ii) a b c d])

    val = 10000*a + 1000*b + 100*c + 10*d + e;
    v = [v; val a b c d e];
    if val==14440
        stim = [stim;onsets(ii) 2];
    elseif val == 12240
        stim = [stim; onsets(ii) 3];
    elseif val == 14400
        stim = [stim; onsets(ii) 4];
    elseif val == 12200 
        stim = [stim; onsets(ii) 5];
    elseif val == 14240
        stim = [stim; onsets(ii) 6];
    elseif val == 12040
        stim = [stim; onsets(ii) 7];
    elseif val == 14200
        stim = [stim; onsets(ii) 8];
    elseif val == 12000
        stim = [stim; onsets(ii) 9];
    elseif val == 14040
        stim = [stim; onsets(ii) 10];
    elseif val == 10440
        stim = [stim; onsets(ii) 11];
    elseif val == 14000
        stim = [stim; onsets(ii) 12];
    elseif val == 10400
        stim = [stim; onsets(ii) 13];
    elseif val == 10000
        stim = [stim; onsets(ii) 1];
    else
        disp('I have no idea what that code is for!!!!')
    end
    
end

if ~isempty(stim)
figure(99)
hist(stim(:,2),-1000:150); title(['stim numbers sorted  number of trials:',num2str(length(stim)/9)])
figure(101)
plot(x)
hold on
plot(stim(:,1),x(stim(:,1)),'*r')
else
    stim=zeros(2,2);
end

function [stim] = getMTplaidtypeIIstim();

global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

z = x;
x = diff(x);

x = resamplevec(x,samp);
z = resamplevec(z,samp);
x(x<0) = 0;
z(z<0) = 0;
x = x./max(x);
z = z./max(z);
z = round(z * 5);
x = round(x * 5);
nx = length(x);
rawpeak = find((z(2:nx-1)==1)&(z(1:nx-2) == 0) & z(3:nx) == 1);
%drop things close to 0
%rawpeak = rawpeak(x(rawpeak)>.7);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
v = [];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 
    
    
    % Let's read out the code for A, B, C, D
    a = 1;
    b= mean(z(onsets(ii)+(8)));
    c = mean(z(onsets(ii)+(14:15)));
    d = mean(z(onsets(ii)+(20:21)));
    e = mean(z(onsets(ii)+(40:41)));

%    disp([onsets(ii) a b c d])

    val = 10000*a + 1000*b + 100*c + 10*d + e;
    v = [v; val a b c d e];
    if val==14440
        stim = [stim;onsets(ii) 2];
    elseif val == 12240
        stim = [stim; onsets(ii) 3];
    elseif val == 14400
        stim = [stim; onsets(ii) 4];
    elseif val == 12200 
        stim = [stim; onsets(ii) 5];
    elseif val == 14240
        stim = [stim; onsets(ii) 6];
    elseif val == 12040
        stim = [stim; onsets(ii) 7];
    elseif val == 14200
        stim = [stim; onsets(ii) 8];
    elseif val == 12000
        stim = [stim; onsets(ii) 9];
    elseif val == 14040
        stim = [stim; onsets(ii) 10];
    elseif val == 10440
        stim = [stim; onsets(ii) 11];
    elseif val == 14000
        stim = [stim; onsets(ii) 12];
    elseif val == 10400
        stim = [stim; onsets(ii) 13];
    elseif val == 10000
        stim = [stim; onsets(ii) 1];
    else
        disp('I have no idea what that code is for!!!!')
    end
    
end

if ~isempty(stim)
figure(99)
hist(stim(:,2),-1000:150); title(['stim numbers sorted  number of trials:',num2str(length(stim)/9)])
figure(101)
plot(x)
hold on
plot(stim(:,1),x(stim(:,1)),'*r')
else
    stim=zeros(2,2);
end

function [stim] = getMTSFtuning();

global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

z = x;
x = diff(x);

x = resamplevec(x,samp);
z = resamplevec(z,samp);
x(x<0) = 0;
z(z<0) = 0;
x = x./max(x);
z = z./max(z);
z = round(z * 5);
x = round(x * 5);
nx = length(x);
rawpeak = find((z(2:nx-1)==1)&(z(1:nx-2) == 0) & z(3:nx) == 1);
%drop things close to 0
%rawpeak = rawpeak(x(rawpeak)>.7);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
v = [];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 
    
    
    % Let's read out the code for A, B, C, D
    a = 1;
    b= mean(z(onsets(ii)+(8)));
    c = mean(z(onsets(ii)+(14:15)));
    d = mean(z(onsets(ii)+(20:21)));
    e = mean(z(onsets(ii)+(40:41)));

%    disp([onsets(ii) a b c d])

    val = 10000*a + 1000*b + 100*c + 10*d + e;
    v = [v; val a b c d e];
    if val==14440
        stim = [stim;onsets(ii) 2];
    elseif val == 12240
        stim = [stim; onsets(ii) 3];
    elseif val == 14400
        stim = [stim; onsets(ii) 4];
    elseif val == 12200 
        stim = [stim; onsets(ii) 5];
    elseif val == 14240
        stim = [stim; onsets(ii) 6];
    elseif val == 12040
        stim = [stim; onsets(ii) 7];
    elseif val == 14200
        stim = [stim; onsets(ii) 8];
    elseif val == 12000
        stim = [stim; onsets(ii) 9];
    elseif val == 14040
        stim = [stim; onsets(ii) 10];
    elseif val == 10440
        stim = [stim; onsets(ii) 11];
    elseif val == 10000
        stim = [stim; onsets(ii) 1];
    else
        disp('I have no idea what that code is for!!!!')
    end
    
end

if ~isempty(stim)
figure(99)
hist(stim(:,2),-1000:150); title(['stim numbers sorted  number of trials:',num2str(length(stim)/9)])
figure(101)
plot(x)
hold on
plot(stim(:,1),x(stim(:,1)),'*r')
else
    stim=zeros(2,2);
end

function [stim] = getMTTFtuning();

global vdat
global green

x = cell2mat(vdat);
samp = size(green,2)/length(x);

z = x;
x = diff(x);

x = resamplevec(x,samp);
z = resamplevec(z,samp);
x(x<0) = 0;
z(z<0) = 0;
x = x./max(x);
z = z./max(z);
z = round(z * 5);
x = round(x * 5);
nx = length(x);
rawpeak = find((z(2:nx-1)==1)&(z(1:nx-2) == 0) & z(3:nx) == 1);
%drop things close to 0
%rawpeak = rawpeak(x(rawpeak)>.7);
%stimulus onsets
onsets=rawpeak(1);
for ii=2:length(rawpeak)
   if (rawpeak(ii)-rawpeak(ii-1))>=40
    onsets = [onsets; rawpeak(ii)];
   end
end

stim=[];
v = [];
ipi = 12; %actually closer to 10 frames
for ii=1:length(onsets) 
    
    
    % Let's read out the code for A, B, C, D
    a = 1;
    b= mean(z(onsets(ii)+(8)));
    c = mean(z(onsets(ii)+(14:15)));
    d = mean(z(onsets(ii)+(20:21)));
    e = mean(z(onsets(ii)+(40:41)));

%    disp([onsets(ii) a b c d])

    val = 10000*a + 1000*b + 100*c + 10*d + e;
    v = [v; val a b c d e];
    if val==14440
        stim = [stim;onsets(ii) 2];
    elseif val == 12240
        stim = [stim; onsets(ii) 3];
    elseif val == 14400
        stim = [stim; onsets(ii) 4];
    elseif val == 12200 
        stim = [stim; onsets(ii) 5];
    elseif val == 14240
        stim = [stim; onsets(ii) 6];
    elseif val == 12040
        stim = [stim; onsets(ii) 7];
    elseif val == 14200
        stim = [stim; onsets(ii) 8];
    elseif val == 12000
        stim = [stim; onsets(ii) 9];
    elseif val == 14040
        stim = [stim; onsets(ii) 10];
    elseif val == 10000
        stim = [stim; onsets(ii) 1];
    else
        disp('I have no idea what that code is for!!!!')
    end
    
end

if ~isempty(stim)
figure(99)
hist(stim(:,2),-1000:150); title(['stim numbers sorted  number of trials:',num2str(length(stim)/9)])
figure(101)
plot(x)
hold on
plot(stim(:,1),x(stim(:,1)),'*r')
else
    stim=zeros(2,2);
end


function [stim] = getMarmieStim();
global vdat
global green

x = cell2mat(vdat);

z = x;
x = diff(x);

x(x<0) = 0;
x =x./(max(x));

onsetsraw = find(x ==1);
samp = length(x)/size(green,2);
onsets = ceil(onsetsraw./samp);
stim(1:length(onsets),1) = onsets;
stim(1:length(onsets),2) = 1;
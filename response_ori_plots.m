

%inds0 = find(data.stim(:,2)==0);
inds1 = find(data.stim(:,2)==1);
inds2 = find(data.stim(:,2)==2);

inds3 = find(data.stim(:,2)==3);
inds4 = find(data.stim(:,2)==4);
inds5 = find(data.stim(:,2)==5);
inds6 = find(data.stim(:,2)==6);
inds7 = find(data.stim(:,2)==7);
inds8 = find(data.stim(:,2)==8);
inds9 = find(data.stim(:,2)==9);
% inds10 = find(data.stim(:,2)==10);

% times0 = data.stim(inds0,1);
times1 = data.stim(inds1,1);
times2 = data.stim(inds2,1);
times3 = data.stim(inds3,1);
times4 = data.stim(inds4,1);
times5 = data.stim(inds5,1);
times6 = data.stim(inds6,1);
times7 = data.stim(inds7,1);
times8 = data.stim(inds8,1);
times9 = data.stim(inds9,1);
% times10 = data.stim(inds10,1);

for b=[2 36 37] %cell number
    w=50; % blank dark cell
    n=5; %repeat number
    
for p=1:n %repeat number
    a = 1;
% f0(p,:) = data.Fraw(b,times0(p)+120:times0(p)+340)-data.Fraw(w,times0(p)+120:times0(p)+340);  
%f0(p,:) = data.neuropil(b,times0(p)+120:times0(p)+340)-data.Fraw(w,times0(p)+120:times0(p)+340);  
f1(p,:) = data.Fraw(b,times1(p)+100:times1(p)+300)-data.Fraw(w,times1(p)+100:times1(p)+300); 
f2(p,:) = data.Fraw(b,times2(p)+100:times2(p)+300)-data.Fraw(w,times2(p)+100:times2(p)+300);
f3(p,:) = data.Fraw(b,times3(p)+100:times3(p)+300)-data.Fraw(w,times3(p)+100:times3(p)+300);
f4(p,:) = data.Fraw(b,times4(p)+100:times4(p)+300)-data.Fraw(w,times4(p)+100:times4(p)+300);
f5(p,:) = data.Fraw(b,times5(p)+100:times5(p)+300)-data.Fraw(w,times5(p)+100:times5(p)+300);
f6(p,:) = data.Fraw(b,times6(p)+100:times6(p)+300)-data.Fraw(w,times6(p)+100:times6(p)+300);
f7(p,:) = data.Fraw(b,times7(p)+100:times7(p)+300)-data.Fraw(w,times7(p)+100:times7(p)+300);
f8(p,:) = data.Fraw(b,times8(p)+100:times8(p)+300)-data.Fraw(w,times8(p)+100:times8(p)+300);
f9(p,:) = data.Fraw(b,times9(p)+100:times9(p)+300)-data.Fraw(w,times9(p)+100:times9(p)+300);
% f10(p,:) = data.Fraw(b,times10(p)+10:times10(p)+170);
end


% for p=1:10
% comb0(p,:)=f0(p,:);
% end

% delta F 
% for p=1:n
%  comb0(p,:)=(f0(p,:)-mean(f0(p,100:115)))./mean(f0(p,100:115));
%  end
 for p=1:n
 comb1(p,:)=(f1(p,:)-mean(f1(p,1:115)))./mean(f1(p,1:115));
 end
for p=1:n
comb2(p,:)=(f2(p,:)-mean(f2(p,1:115)))./mean(f2(p,1:115));
end
for p=1:n
comb3(p,:)=(f3(p,:)-mean(f3(p,1:115)))./mean(f3(p,1:115));
end
for p=1:n
comb4(p,:)=(f4(p,:)-mean(f4(p,1:115)))./mean(f4(p,1:115));
end
for p=1:n
comb5(p,:)=(f5(p,:)-mean(f5(p,1:115)))./mean(f5(p,1:115));
end
for p=1:n
comb6(p,:)=(f6(p,:)-mean(f6(p,1:115)))./mean(f6(p,1:115));
end
for p=1:n
comb7(p,:)=(f7(p,:)-mean(f7(p,1:115)))./mean(f7(p,1:115));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:n
comb8(p,:)=(f8(p,:)-mean(f8(p,1:115)))./mean(f8(p,1:115));
end

for p=1:n
comb9(p,:)=(f9(p,:)-mean(f9(p,1:115)))./mean(f9(p,1:115));
end

%standard dev
% for p=1:n
% comb1(p,:)=(f1(p,:)-mean(mean(f1)))./mean(std(f1(p,:)));
% end
%  for p=1:n
% comb2(p,:)=(f2(p,:)-mean(mean(f1(p,:))))./mean(std(f1(p,:)));
%  end
% for p=1:n
% comb3(p,:)=(f3(p,:)-mean(mean(f1(p,:))))./mean(std(f1(p,:)));
% end
% for p=1:n
% comb4(p,:)=(f4(p,:)-mean(mean(f1(p,:))))./mean(std(f1(p,:)));
% end
% for p=1:n
% comb5(p,:)=(f5(p,:)-mean(mean(f1(p,:))))./mean(std(f1(p,:)));
% end
% for p=1:n
% comb6(p,:)=(f6(p,:)-mean(mean(f1(p,:))))./mean(std(f1(p,:)));
% end
% for p=1:n
% comb7(p,:)=(f7(p,:)-mean(mean(f1(p,:))))./mean(std(f1(p,:)));
% end







%comb0 = comb0;

comb = [comb1,comb2,comb3,comb4,comb5,comb6,comb7,comb8,comb9];
combFraw = [f1,f2,f3,f4,f5,f6,f7,f8,f9];
combFraw = reshape(combFraw,5,201,9);

combf = reshape(comb,5,201,9); %5 repeat, 141 frame 9 stim


combfa = mean(combf,1);
combFrawa = mean(combFraw,1);
%figure(b+100)
%for j=1:9
%subplot(3,3,j)        
%plot(combfa(:,:,j))

   %end
 
time(1) = 1/32;
for i =2: 201
    time(i) = time(i-1) + 1/32;
end
   
figure (b+300)

for j=1:9 %num of stim. condition
subplot(1,9,j)
   for i=1:n          
plot(time,combf(i,:,j),'b')
hold on
plot(time,combfa(:,:,j),'r')
ylim([-2 2])
xlim([0 6])
line([0.625 0.625],[0 2]) %stim onset
hold on
%line([4 4],[0 1])

   end
   
end

% figure (b+300)
% 
% for j=1:9 %num of stim. condition
% subplot(9,1,j)
%    for i=1:n          
% plot(time,combFraw(i,:,j),'k')
% hold on
% plot(time,combFrawa(:,:,j),'r')
% %ylim([-1 1])
% %xlim([0 6])
% %line([1 1],[0 1]) %stim onset
% hold on
% %line([3 3],[0 1])
% 
%    end
% %    
%  end


end
function varargout = Caliman(varargin)
% GUI to reads calcium imaging files from labview
% Sari Andoni
% Priebe Lab
% University of Texas

% Last Modified by GUIDE v2.5 03-Sep-2016 19:36:41

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Caliman_OpeningFcn, ...
                   'gui_OutputFcn',  @Caliman_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Caliman is made visible.
function Caliman_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Caliman (see VARARGIN)

% Choose default command line output for Caliman
handles.output = hObject;
vv = ver('MATLAB');
%Set callback for continuous slider
if str2num(vv.Version)>8.5
hListener=addlistener(handles.slider1,'ContinuousValueChange',@slider1_Callback);
else

hListener=handle.listener(handles.slider1,'ActionEvent',@slider1_Callback);
end

%setappdata(handles.slider1,'myListener',hListener);
handles.hListener=hListener;

handles.oldschool = 0;

handles.isPlaying=false;
handles.avgOn=0;
handles.neuronIndex=0;
handles.mask=[];
handles.addMode=1;
handles.colors=[0         0    1.0000;
                1    .5         0.5;
            1.0000         0         0;
                 0    0.7500    0.7500;
            0.7500         0    0.7500;
            0.7500    0.7500         0;
            0.2500    0.7500    0.2500];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Caliman wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Caliman_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_ButtonDownFcn(hObject, eventdata, handles)
return


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles=guidata(hObject);
handles.avgOn=0;
x=handles.nx;
y=handles.ny;
%handles
f=floor(get(hObject,'Value'));
updateimage(f,handles);
% v=get(handles.popupmenu1,'Value');
% global green
% global red
% if v==1
%     raw=squeeze(green{f});
%     raw=raw-min(raw(:));
%     raw=raw/max(raw(:));
%     img=cat(3,raw,raw,raw);
% elseif v==2
%     raw=squeeze(red{f});
%     raw=raw-min(raw(:));
%     raw=raw/max(raw(:));
%     img=cat(3,raw,raw,raw);
% else
%     raw1=squeeze(green{f});
%     raw2=squeeze(red{f});
%     raw2=raw2-min(raw2(:));
%     raw2=raw2/max(raw2(:));
%     raw1=raw1-min(raw1(:));
%     raw1=raw1/max(raw1(:));
%     img=cat(3,raw2,raw1,zeros(x,y));
% end
%     
% set(handles.rimg,'CData',img);
dframe = round(get(handles.slider1,'Value'));
set(handles.framenumber,'String',dframe);
if strcmp(handles.file(20:21),'DS')
    % list depth
    de = round(1000*((dframe/handles.length) * (handles.desc.Ztop - handles.desc.Zbot) + handles.desc.Zbot));
    set(handles.time,'String',de);
else 
    set(handles.time,'String',round(100*round(get(handles.slider1,'Value'))*256/7910)/100);
end
title(handles.img,strcat('Frame #',num2str(f)));
set(handles.rimg, 'ButtonDownFcn', {@img_ButtonDownFcn, handles});
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% 
% % --- Executes on slider movement.
% function slider2_Callback(hObject, eventdata, handles)
% % hObject    handle to slider2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'Value') returns position of slider
% %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% plot(handles.img,randn(100,1));
% drawnow;
% 
% % --- Executes during object creation, after setting all properties.
% function slider2_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to slider2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: slider controls usually have a light gray background.
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end


% --- Executes on button press in findCellsBtn.
function findCellsBtn_Callback(hObject, eventdata, handles)
% hObject    handle to findCellsBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sprintf('Finding Cells....');

% --- Executes on button press in regularizeBtn.
function regularizeBtn_Callback(hObject, eventdata, handles)
% hObject    handle to regularizeBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles.rgreen=handles.green;
axes(handles.img);
%cla;
a=axis;
li = text((a(2)-a(1))/2,(a(4)-a(3))/2,'Regularizing...');
set(li,'Color',[0 0 1])
drawnow;
zout = fliplr([45 47 124 92 ]);
zoutind = 0;

global green red
 %img1=double(squeeze(handles.green(:,:,1)));
% img1=handles.gmImg;    
 img1=double(green{1});
 mrimg = double(red{1});
%handles
 for i=2:handles.length,
     img2=double(green{i});
     [p,newImg]=dftregistration(fft2(img1),fft2(img2),10);
     green{i}=real(ifft2(dftshift(fft2(green{i}),p(2),p(3),p(end))));
     red{i}=real(ifft2(dftshift(fft2(red{i}),p(2),p(3),p(end))));
%     green{i}=real(ifft2(newImg));
     img1 = img1+double(green{i});

     mrimg = mrimg + double(red{i});
     handles.offsets(i,:) = handles.offsets(i,:) + p(3:4);
     fprintf('%c',[8 zout(rem(i,4)+1)])
     zoutind = 1-zoutind;
 end
 
 handles.gmImg2 = img1/handles.length;
 
 handles.rmImg2 = mrimg/handles.length;
 delete(li)
 updateimage(floor(get(handles.slider1,'Value')),handles);
 
 disp('Regularized!')
% handles.rimg=imagesc(handles.mImg,'Parent',handles.img);
% set(handles.rimg, 'ButtonDownFcn', {@img_ButtonDownFcn, handles});
% axis(handles.img,'xy');
guidata(hObject,handles)



% --- Executes on button press in openFile.
function openFile_Callback(hObject, eventdata, handles)
% hObject    handle to openFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if exist('lastFile.mat','file'),
    load('lastFile.mat');
else
    lastfile='';
end
handles.output = hObject;
vv = ver('MATLAB');
if str2num(vv.Version)>8.5
[file,path]=uigetfile('*.mat','',lastfile);
else
[file,path,filterindex]=uigetfile({'*.mat'},'',lastfile);
end
if length(file)==1 && file==0
    return
end

lastfile=strcat(path,file);
% save('lastFile.mat','lastfile');
axes(handles.img);
cla;
a=axis;
text((a(2)-a(1))/2,(a(4)-a(3))/2,'Loading...');
drawnow;

numFrames=0;%str2num(get(handles.numFrames,'String'));

if numFrames > 0,
%     l2phfile_fast(handles,path,file,0,numFrames);
%     l2phfile_resFF(handles,path,file,0,numFrames);
    [desc, gccode] = l2phfile_resFF2(path,file,0,1:numFrames);
else
%     l2phfile_fast(handles,path,file);
%     l2phfile_resFF(handles,path,file);
%    [desc, gccode] = l2phfile_resFF2(path,file);
     load([lastfile]);
end
global green
global red
global vdat
if ~isfield(desc(1),'Zoom')
    desc(1).Zoom = [300 300];  % Let this be default for "old" data
else
    desc(1).Zoom = [desc(1).Zoom desc(1).Zoom];
end

if handles.oldschool
    disp('OLDSCHOOL!')
    desc(1).Zoom = [455 256];
end

handles.desc = desc;

% Check to see if it is bidirectional:
datasize = size(red{1});

if datasize(2)>500  %% This is bidirectional
    prompt='This is bidirectional data.  Do you want  dual or single direction?';
    operation = 'Bidirectional handling';
    choice2 = 'Dual';
    choice1 = 'Single';
    choice = questdlg(prompt, operation, choice1, choice2, choice1);
    hz = datasize(2)/2;
    if strcmp(choice,choice2)
      prompt='Do you want  interleave or align right and left signals?';

      operation = 'Bidirectional integration';
      choice1 = 'Interleave';
      choice2 = 'Align';
      choiceint = questdlg(prompt, operation, choice1, choice2, choice1);

      if strcmp(choiceint,choice1)  %% No alignment, just flip images
        a = zeros(size(red{1},1)*2,size(red{1},2)/2);
        for j=1:length(red)
          a(1:2:end) = red{j}(:,1:hz);
          a(2:2:end) = fliplr(red{j}(:,hz+1:end));
          red{j} = a;
          a(1:2:end) = green{j}(:,1:hz);
          a(2:2:end) = fliplr(green{j}(:,hz+1:end));
          green{j} = a;
        end
      else  % Align and flip
        prompt='Which channel do you want to use for alignment?';
        operation = 'Channel';
        choice1 = 'red';
        choice2 = 'green';
        choicech = questdlg(prompt, operation, choice1, choice2, choice1);
        

        step = 100;  %% Could be varied...
        cnt =1;
        lstep = floor(size(red,2)/step)*step;
        rotamnt = [];
        for j=1:step:lstep
            rr = zeros(size(red{1},1),hz);
            rr1 = zeros(size(rr));
            if strcmp(choicech,choice1)
              for k=1:step
                rr= rr+red{(j-1)+k}(:,1:hz);
                rr1 = rr1+fliplr(red{(j-1)+k}(:,hz+1:end));
              end
            else
              for k=1:step
                rr= rr+green{(j-1)+k}(:,1:hz);
                rr1 = rr1+fliplr(green{(j-1)+k}(:,hz+1:end));
              end
            end
            % Now we have some data, let's align.
            [p]=dftregistration(fft2(rr),fft2(rr1),10);
            rotamnt(cnt) = p(end);
            dph(cnt) = p(2);
            cnt = cnt + 1;
        end
        rotamnt = median(rotamnt);
        dph = median(dph);
        %% OK, now we have the rotation step size.  We need to do it now for all the frames:
        fprintf('Alignment offset = %d\n',dph);
        a = zeros(size(red{1},1)*2,size(red{1},2)/2);
        for j=1:size(red,2)
            a(1:2:end) = red{j}(:,1:hz);
            a(2:2:end) = real(ifft2(dftshift(fft2(fliplr(red{j}(:,hz+1:end))),dph,0,rotamnt)));
            red{j} = a;
            
            a(1:2:end) = green{j}(:,1:hz);
            a(2:2:end) = real(ifft2(dftshift(fft2(fliplr(green{j}(:,hz+1:end))),dph,0,rotamnt)));
            green{j} = a;

        end
        
      end
    else
        operation = 'Direction';
        choice1 = 'left';
        choice2 = 'right';
        choicedir = questdlg(prompt, operation, choice1, choice2, choice1);
        if strcmp(choicedir, choice1)
          for j=1:length(red)
            red{j} = red{j}(:,1:hz);% + fliplr(red{j}(:,hz+1:end));
            green{j} = green{j}(:,1:hz);% + fliplr(green{j}(:,hz+1:end));
          end
        else
          for j=1:length(red)
            red{j} = fliplr(red{j}(:,hz+1:end));
            green{j} = fliplr(green{j}(:,hz+1:end));
          end
        end
    end
 
end    
handles.file=file;
handles.offsets = zeros(size(green,2),2);


% %%%%%%%%%%%%%%%
% % STEADY CAM2
% %params to update for each imaging session
% %how to identify bright thing automatically?
% 
% %%%
% thresh = .5;
% win1 = 120:170;
% win2 = 200:250;
% %%%
% 
% gwid = 1;
% gaus = exp(-((-2:2).^2)./(2*gwid*gwid));
% gaus2d = gaus'*gaus;
% 
% template = green{500};
% [origs1,origs2] = size(template);
% template = template(win1,win2);
% template = template./max(max(template));
% template(template<thresh)=0;
% template = conv2(template,gaus2d,'same');
% [s1,s2] = size(template);
% for fr=1:5900
%     %tic
%     zerostemplate = zeros(origs1+40,origs2+80); 
%     %add 10 pixels and 20 pixels to locations, then do xcorr adjustment
%     offsetTemplate = green{fr};
%     offsetTemplate = offsetTemplate(win1,win2);
%     offsetTemplate = offsetTemplate./max(max(offsetTemplate));
%     offsetTemplate(offsetTemplate<thresh)=0;
%     offsetTemplate = conv2(offsetTemplate,gaus2d,'same');
%     cc = xcorr2(offsetTemplate,template);
%     [max_cc, imax] = max(abs(cc(:)));
%     [ypeak, xpeak] = ind2sub(size(cc),imax(1));
%    
%     xoffset = ypeak-s1;
%     yoffset = xpeak-s2;
%     
%     xl = (20:origs1+20-1) - xoffset;
%     yl = (40:origs2+40-1) - yoffset;
%     
%     zerostemplate(xl,yl) = green{fr};
%     green{fr} = zerostemplate;
%     fr
% end
% %%%%%%%%%%%%%
 
% stim = getORIstim_galvos(data);
%stim = getORIstim(); 
%%%%%%%%%%%%%%%%%%%%%%%%stim = getORIstim2(); 


if strcmp(file(20:21),'CC')  %% Single plane acquisition
    [handles.stim, handles.gcstim handles.frcode] = getstim(gccode);   %% NEEDS TO BE MORE GENERAL!
else
    handles.stim = [];
end
% if strcmp(file(20:21),'CC')  %% Single plane acquisition
%     [handles.stim, handles.gcstim] = getstim_old(gccode);   %% NEEDS TO BE MORE GENERAL!
% else
%     handles.stim = [];
% end
    
%stim = getBINOCstim2();

% stim(:,1)=ceil(stim(:,1));
% data.stim = stim;


%number frames to avg
% s =[];
% for ii = 1:size(stim,1)
%     s = [s;  stim(ii,1)+ 100:stim(ii,1)+ 200];
% end
% nnn = [];
% for ii = 1:size(stim,1)*size(stim,2)
%     nnn = [nnn; stim(ii)];
% end
% 
% 
%  s = [];
nnn=length(green);
g1=double(green{1});
[x,y]=size(g1);
handles.nx=x;
handles.ny=y;
l=length(green);
handles.gmImg=g1;
handles.gmImg2=g1;
for ii=2:nnn
handles.gmImg = handles.gmImg+double(green{ii});
end
handles.gmImg = handles.gmImg./nnn;
handles.gmImg2 = handles.gmImg;
handles.gmImg=handles.gmImg - min(handles.gmImg(:));
%handles.gmImg2 = handles.gmImg;
handles.gmImg=handles.gmImg / max(handles.gmImg(:));
handles.length=l;
handles.mImg=cat(3,handles.gmImg,handles.gmImg,handles.gmImg);

r1=double(red{1});
handles.rmImg=r1;
handles.rmImg2=r1;
for ii=2:nnn
handles.rmImg = handles.rmImg+double(red{ii});
end
handles.rmImg = handles.rmImg ./nnn;
handles.rmImg2 = handles.rmImg;
handles.rmImg=handles.rmImg - min(handles.rmImg(:));
%handles.rmImg2 = handles.rmImg;
handles.rmImg=handles.rmImg / max(handles.rmImg(:));


xs = desc.Zoom(1) *(1:handles.ny)./handles.ny;
ys = desc.Zoom(2) *(1:handles.nx)./handles.nx;
axes(handles.img);
handles.rimg=imagesc(xs,ys,handles.mImg);
set(handles.img,'Box','off');
set(handles.meanimage,'Value',1);
axis xy;
handles.avgOn=1;


set(handles.slider1,'Value',1);
set(handles.slider1,'Max',handles.length);
set(handles.slider1,'Min',1);
set(handles.text1,'String',lastfile);
set(handles.rimg, 'ButtonDownFcn', {@img_ButtonDownFcn, handles});
set(handles.popupmenu1,'Value',1);

set(handles.framenumber,'String',round(get(handles.slider1,'Value')));
set(handles.time,'String',round(100*round(get(handles.slider1,'Value'))*256/7910)/100);
handles.neuron=[];
handles.mask=[];
handles.neuronIndex=0;

% Update handles structure
guidata(hObject, handles);






% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cind = floor(get(handles.slider1,'Value'));
while (get(handles.play,'Value') == 1) & (cind <  (get(handles.slider1,'Max')))
    cind = 1+cind;
    v=get(handles.popupmenu1,'Value');
    x=handles.nx;
    y=handles.ny;
    global green
    global red
    updateimage(cind,handles)
    set(handles.slider1,'Value',cind);
    drawnow;
    set(handles.framenumber,'String',round(get(handles.slider1,'Value')));
    set(handles.time,'String',round(100*round(get(handles.slider1,'Value'))*256/7910)/100);

    set(handles.rimg, 'ButtonDownFcn', {@img_ButtonDownFcn, handles});
    pause(0.00000005)
    
end


% --- Executes on button press in average.
function average_Callback(hObject, eventdata, handles)
% hObject    handle to average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles.rimg=imagesc(handles.mImg,'Parent',handles.img);
%axis(handles.img,'xy');
set(handles.rimg,'CData',handles.mImg);
set(handles.rimg, 'ButtonDownFcn', {@img_ButtonDownFcn, handles});
handles.avgOn=1;
%    handles.green=handles.rgreen;
guidata(hObject, handles);

% --- Executes on button press in roi.
function roi_Callback(hObject, eventdata, handles)
% hObject    handle to roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi


% --- Executes on button press in pixel.
function pixel_Callback(hObject, eventdata, handles)
% hObject    handle to pixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pixel


% --- Executes on button press in cell.
function cell_Callback(hObject, eventdata, handles)
% hObject    handle to cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell

% --- Executes on button press in cell.
function image_Callback(hObject, eventdata, handles)
% hObject    handle to cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cind = floor(get(handles.slider1,'Value'));

updateimage(cind, handles);


% --- Executes on mouse press over axes background.
function img_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=get(hObject,'Parent');
handles=guidata(h);
handles.rcont=[];
pos = get(handles.img,'Currentpoint');
pos =(pos(1,1:2));

%handles.rimg=imagesc(handles.mImg,'Parent',handles.img);
axis(handles.img,'xy');
handles.addMode

if handles.addMode == 1,

    pos(1) = round(pos(1)*handles.ny/handles.desc.Zoom(1));
    pos(2) = round(pos(2)*handles.nx/handles.desc.Zoom(2));
    W2blck = str2double(get(handles.cellScale,'String'));
    W2bh = floor(W2blck/2);

    [x,y] = meshgrid(-W2bh:W2bh,-W2bh:W2bh);
    r = sqrt(x.^2 + y.^2);
    sig = W2blck/3;
    Win = exp(-r.^2/(2*sig^2));

    xx=pos(2)-W2bh:pos(2)+W2bh;
    yy=pos(1)-W2bh:pos(1)+W2bh;
    
    if sum(xx<1)>0
        Win=Win(xx>0,xx>0);
        yy=yy(xx>0);
        xx=xx(xx>0);
    elseif sum(yy<1)>0
        Win=Win(yy>0,yy>0);
        xx=xx(yy>0);
        yy=yy(yy>0);
    end
    
    if xx(end)>=handles.nx
        Win=Win(xx<=handles.nx,xx<=handles.nx);
        yy=yy(xx<=handles.nx);
        xx=xx(xx<=handles.nx);
    end
    
    if yy(end)>=handles.ny
        Win=Win(yy<=handles.ny,yy<=handles.ny);
        xx=xx(yy<=handles.ny);
        yy=yy(yy<=handles.ny);
    end
    
    
    blck = handles.gmImg2(xx,yy);
    blck = blck-min(blck(:));
    blck = blck/max(blck(:));
    blck = blck.*Win;
    idcell = find(blck>.3);

    bwCell=zeros(size(handles.gmImg));

    blckBW = zeros(size(blck)); 
    blckBW(idcell) = 1;
    blck = bwCell(xx,yy);
    blck = sign(blck + blckBW);
    bwCell(xx,yy) = blck;

    handles.neuronIndex=handles.neuronIndex+1;
    handles.mask(:,:,handles.neuronIndex)=bwCell;
    
elseif handles.addMode == 2,
    pos(1) = round(pos(1)*handles.ny/handles.desc.Zoom(1));
    pos(2) = round(pos(2)*handles.nx/handles.desc.Zoom(2));
%pos %= round(pos)
    removed=false;
%    fprintf('Hey: %g\n',handles.neuronIndex)
    for i=1:handles.neuronIndex,
        if handles.mask(pos(2),pos(1),i)>0,
            handles.mask=cat(3,handles.mask(:,:,1:i-1),handles.mask(:,:,i+1:end));
            handles.neuronIndex=handles.neuronIndex-1;
            removed=true;
%            disp('removed')
            break;
        end
    end
    
    if ~removed,
        beep;
    end
    
else
    axes(handles.img);
    bwCell=roipoly;
    %bwCell=flipud(bwCell');
    handles.neuronIndex=handles.neuronIndex+1;
    handles.mask(:,:,handles.neuronIndex)=bwCell;
end
imgch = get(handles.img,'Children');
delete(imgch(1:end-1));

if handles.neuronIndex > 0,

hold(handles.img, 'on');
l=handles.length;
xs = handles.desc.Zoom(1) *(1:handles.ny)./handles.ny;
ys = handles.desc.Zoom(2) *(1:handles.nx)./handles.nx;
for i=1:handles.neuronIndex,
    c=handles.colors(mod(i-1,7)+1,:);
    contour(handles.img,xs,ys,handles.mask(:,:,i),2,'Color',c,'HitTest','off');
    hold on
    [y,x]=find(handles.mask(:,:,i) > 0);
    y=median(y)*handles.desc.Zoom(2)./handles.nx;x=median(x)*handles.desc.Zoom(1)./handles.ny;
    text(x,y,num2str(i),'Color',c,'FontSize',10,'Parent',handles.img);
end
hold(handles.img, 'off');

end

set(handles.rimg, 'ButtonDownFcn', {@img_ButtonDownFcn, handles});

guidata(h,handles);






% --- Executes on button press in saveAs.
function saveAs_Callback(hObject, eventdata, handles)
% hObject    handle to saveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles;
prompt='Do you want to write the Red or Green channel?';
operation = 'Channel';
choice2 = 'Green';
channelchoice = 'Red';
channelchoice = questdlg(prompt, operation, channelchoice, choice2, channelchoice);

[file,path]=uiputfile('*.tiff');
sx=strfind(file,'.');
ext=file(sx+1:end);
prefix = file(1:(sx-1));

prompt={'Enter the number of frames to average (1=none)'};
operation='Bins';
numlines= 1;
defaultanswer={'1'};
answer=inputdlg(prompt,operation,numlines,defaultanswer);

prompt='Do you want the CLUT to be fixed or variable?';
operation = 'CLUT';
clutchoice = 'Fixed';
choice2 = 'Variable';
clutchoice = questdlg(prompt, operation, clutchoice, choice2, clutchoice);



framestep = max([str2num(answer{1}) 1]);
firsttime = 1;
if strcmp(channelchoice, 'Red') % Red
    global red
    
    if (size(red{1},1)==512)
        ym = 500;
    else
        ym = 250;
    end
    
    minv = round(min(red{end}(:))); 
    maxv = round(max(red{end}(:)));
    if strcmp(clutchoice,'Fixed')
        prompt={'Min Value'; 'Max Value'};
        operation='CLUT range';
        numlines= 1;
        defaultanswer={num2str(minv),num2str(maxv)};
        clutanswer=inputdlg(prompt,operation,numlines,defaultanswer);
        if length(clutanswer)>0
            minv = str2num(clutanswer{1});
            maxv = str2num(clutanswer{2});
        end
    end     
    for k = 1:framestep:(framestep*floor(handles.length/framestep))
       out = double(red{k}(1:ym,:));
       for j=(k+1):(k+framestep-1)
           out = out + red{j}(1:ym,:);
       end
       
       out = squeeze(out/framestep);
       if strcmp(clutchoice,'Variable')
           minv = min(out(:)); maxv = max(out(:));
       end
       out = out-minv; out = uint8(255*out./maxv);
       if firsttime
           imwrite(flipud(out),(gray(256)),[path file],'tiff','writemode','overwrite','Compression','none');
       else
           imwrite(flipud(out),(gray(256)),[path file],'tiff','writemode','append','Compression','none');
       end
       firsttime = 0;
    end
else
    global green
    if (size(green{1},1)==512)
        ym = 500;
    else
        ym = 250;
    end
    minv = round(min(green{end}(:)));
    maxv = round(max(green{end}(:)));
    if strcmp(clutchoice,'Fixed')
        prompt={'Min Value'; 'Max Value'};
        operation='CLUT range';
        numlines= 1;
        defaultanswer={num2str(minv),num2str(maxv)};
        clutanswer=inputdlg(prompt,operation,numlines,defaultanswer);
        if length(clutanswer)>0
            minv = str2num(clutanswer{1});
            maxv = str2num(clutanswer{2});
        end
    end     
    for k = 1:framestep:(framestep*floor(handles.length/framestep))
       out = double(green{k}(1:ym,:));
       for j=(k+1):(k+framestep-1)
           out = out + green{j}(1:ym,:);
       end
       
       out = squeeze(out/framestep);
       if strcmp(clutchoice,'Variable')
           minv = min(out(:)); maxv = max(out(:));
       end
       out = out-minv; out = uint8(255*out./maxv);
       if firsttime
           imwrite(flipud(out),(gray(256)),[path file],'tiff','writemode','overwrite','Compression','none');
       else
           imwrite(flipud(out),(gray(256)),[path file],'tiff','writemode','append','Compression','none');
       end
       firsttime = 0;
    end
end

% [file,path]=uiputfile('*.avi;*.mat');
% sx=strfind(file,'.');
% ext=file(sx+1:end);
% if strcmp(ext,'avi'),
%     % Prepare the new file.
%     vidObj = VideoWriter(strcat(path,file));
%     open(vidObj);
% 
%     for k = 1:handles.length, 
%        % Write each frame to the file.
%        green=double(handles.green(:,:,k));
%        green=green-min(green(:));
%        green=green./max(green(:))*255;
%        green=uint8(green');
%        writeVideo(vidObj,flipud(green));
%     end
% 
%     % Close the file.
%     close(vidObj);
% % else
% %     green=handles.green;
% %     save(file,'green');
% end

    
    


% --- Executes on button press in saveImage.
function saveImage_Callback(hObject, eventdata, handles)

% hObject    handle to saveImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k=floor(get(handles.slider1,'Value'));
[file,path]=uiputfile('*.jpg');

if get(handles.meanimage,'Value'),
    igreen=handles.mImg;
    igreen=igreen./max(igreen(:))*255;
    igreen=uint8((igreen));
end
imwrite(igreen,strcat(path,file),'jpg');


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

cind = floor(get(handles.slider1,'Value'));
updateimage(cind, handles);

% v=get(hObject,'Value');
% x=handles.nx;
% y=handles.ny;
% if v == 1,
%     handles.mImg=cat(3,handles.gmImg,handles.gmImg,handles.gmImg);
% elseif v==2,
%     handles.mImg=cat(3,handles.rmImg,handles.rmImg,handles.rmImg);
% else
%     handles.mImg=cat(3,handles.rmImg,handles.gmImg,zeros(x,y));
% end
% 
% handles.rimg=imagesc(handles.mImg,'Parent',handles.img);
% axis(handles.img,'xy');
% set(handles.rimg, 'ButtonDownFcn', {@img_ButtonDownFcn, handles});
% handles.avgOn=1;

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addCell.
function addCell_Callback(hObject, eventdata, handles)
% hObject    handle to addCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.addMode=1;
guidata(hObject, handles);

% --- Executes on button press in removeCell.
function removeCell_Callback(hObject, eventdata, handles)
% hObject    handle to removeCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.addMode=2;
guidata(hObject, handles);


% --- Executes on button press in addRoi.
function addRoi_Callback(hObject, eventdata, handles)
% hObject    handle to addRoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.addMode=3;
guidata(hObject, handles);


% % --- Executes on button press in export.
% function export_Callback(hObject, eventdata, handles)
% % hObject    handle to export (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% export2wsdlg({'Cell Responses','Cell Masks','Stimulus','All Frames'},...
%              {'r','mask','stim','frames'},...
%              {handles.neuron,handles.mask,handles.stim,handles.green},...
%              'Export data to Matlab Workspace',...
%              [true,true,true,false]);

%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)

global vdat
global green
global red

% darknoise = -107.6; %% 0.12 based on 13 Jun 2019 %% -107.6 based on 17 April 2018 by JJP; %% before -160; mean for experiments 04/24/2015 to 04/27/2015
% 
% addpath C:\Users\JagrutiJP\Projects\RawData\2phlmatlab\

mm=zeros(1,handles.length);
Fraw=[];
disp(['frames: ',num2str(handles.length)])
disp(['masks: ',num2str(size(handles.mask))])
disp('getting Fraw: ')
%get cell data
for cc = 1:size(handles.mask,3)
    disp(num2str(cc))
    bwCell = logical(squeeze(handles.mask(:,:,cc)));
    for ii=1:length(mm)
        mm(ii) = mean(green{ii}(bwCell));
    end
    %Fraw(cc,:)=mm-darknoise;
    Fraw(cc,:)=mm;
end
 
file=handles.file;
if exist(['\2pht\analyzed\',file(1:10)])<6
    mkdir(['\2pht\analyzed\',file(1:10)])
end
locfile =['\2pht\analyzed\',file(1:10),'\',file(12:18)];
disp(locfile)

data=[];
data.filename = file;
data.masks = handles.mask;
data.sampleG = handles.gmImg2;
data.sampleR = handles.rmImg2;
%data.green = green;
data.vdat = vdat;
data.stim = handles.stim;
data.gcstim = handles.gcstim;
data.frcode = handles.frcode
data.Fraw = Fraw;

%get average pixel number per mask
maskpixelsize=[];
template = zeros(size(handles.gmImg));
for cc=1:size(data.masks,3)
    maskpixelsize(cc) = sum(sum(squeeze(data.masks(:,:,cc))));
    template = template+squeeze(data.masks(:,:,cc));
end
maskpixelsize = round(mean(maskpixelsize));
%rebuild neuropil with reverse mask template
template = ~(template); 
%remove image edges to get neuropil (in case of negative values)
template(1:50,:) = 0;
template(256-50:256,:) = 0;
template(:,1:85) = 0; 
template(:,455-85:455) = 0;
template = logical(template);
%get neuropil data
for ii=1:length(mm)
    mm(ii)=mean(green{ii}(template)); 
end
%data.neuropil = mm-darknoise;
data.neuropil = mm;




save(locfile,'data', '-v7.3')

export2wsdlg({'data'},{'data'},{data},'Export data to Matlab Workspace',[true]);

         
         
% --- Executes on button press in clearCells.
function clearCells_Callback(hObject, eventdata, handles)
% hObject    handle to clearCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.neuron=[];
handles.neuronIndex=0;
handles.mask=[];
imgch = get(handles.img,'Children');
delete(imgch(1:end-1));
%handles.rimg=imagesc(handles.mImg,'Parent',handles.img);
%axis(handles.img,'xy');
%set(handles.rimg, 'ButtonDownFcn', {@img_ButtonDownFcn, handles});
%handles.avgOn=1;

guidata(hObject, handles);



function cellScale_Callback(hObject, eventdata, handles)
% hObject    handle to cellScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cellScale as text
%        str2double(get(hObject,'String')) returns contents of cellScale as a double


% --- Executes during object creation, after setting all properties.
function cellScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cellScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numFrames_Callback(hObject, eventdata, handles)
% hObject    handle to numFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numFrames as text
%        str2double(get(hObject,'String')) returns contents of numFrames as a double


% --- Executes during object creation, after setting all properties.
function numFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadmasks.
function loadmasks_Callback(hObject, eventdata, handles)
% hObject    handle to loadmasks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%



% fileloc=['/Users/benscholl/Documents/2phl/analyzed/',...
%     handles.file(1:10),'/',handles.file(12:18),'.mat'];
[file,path]=uigetfile('*.mat','','');
fileloc =strcat(path,file);
if exist(fileloc)~=0
   load(fileloc)
   h=get(hObject,'Parent');
   handles=guidata(h);
   %% First let's get rid of any old masks:
   ch = get(handles.img,'Children');
   if length(ch)>1
       delete(ch(1:end-1))
   end
   
   masks = data.masks;
   shiftx = 20;%%pos is left
   shifty =15;%% pos is down
   for cc = 1:size(masks,3)
       for xj = 1:size(masks, 1)
           for yk = 1:size(masks, 2)
               
               if (xj+shifty) > 256
                   ysh = rem((xj+shifty), 256);
               elseif (xj+ shifty) <= 0
                   ysh = 256 + (xj + shifty);
               else
                   ysh = (xj+ shifty) ;
               end
               if (yk+shiftx) > 455
                   xsh = rem((yk+shiftx), 455);
               elseif (yk+shiftx) <= 0
                   xsh = 455 + (yk+shiftx);
               else
                   
                   xsh = yk+shiftx;
               end                  
               newmasks(xj, yk, cc) = masks(ysh, xsh, cc);
           end
       end
       %masks(:,:,cc) = circshift(masks(:,:,cc), 2, 50);    
   end 
   handles.mask = newmasks;   %%%% translate mask 
   handles.neuronIndex = size(data.masks,3);
   %handles.stim = data.stim;
   hold(handles.img, 'on');
       
       xs = handles.desc.Zoom(1) *(1:handles.ny)./handles.ny;
       ys = handles.desc.Zoom(2) *(1:handles.nx)./handles.nx;
   
   
   for i=1:handles.neuronIndex,
       c=handles.colors(mod(i-1,7)+1,:);
       contour(handles.img,xs,ys,handles.mask(:,:,i),1,'Color',c,'HitTest','off');
       [y,x]=find(handles.mask(:,:,i) > 0);
       y=median(y) * handles.desc.Zoom(2)/handles.nx;x=handles.desc.Zoom(1) * median(x)/handles.ny;
       text(x,y,num2str(i),'Color',c,'FontSize',10,'Parent',handles.img);
   end
   hold(handles.img, 'off');
   set(handles.rimg, 'ButtonDownFcn', {@img_ButtonDownFcn, handles});
   guidata(h,handles);
else
    disp('No file saved')
end



function [Greg] = dftshift(buf2ft,diffphase,row_shift, col_shift)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007

% Portions of this code were taken from code written by Ann M. Kowalczyk 
% and James R. Fienup. 
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued 
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458 
% (1990).

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).

% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)

% Outputs
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.

% Compute registered version of buf2ft
if (abs(col_shift) > 0)
    [nr,nc]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
else
    Greg = buf2ft;
end
return

function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
return



function [output Greg] = dftregistration(buf1ft,buf2ft,usfac)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007

% Portions of this code were taken from code written by Ann M. Kowalczyk 
% and James R. Fienup. 
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued 
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458 
% (1990).

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).

% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)

% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.

% Default usfac to 1
if exist('usfac')~=1, usfac=1; end

% Compute error for no pixel shift
if usfac == 0,
    CCmax = sum(sum(buf1ft.*conj(buf2ft))); 
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2); 
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero); 
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax)); 
    output=[error,diffphase];
        
% Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
% peak
elseif usfac == 1,
    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);
    cloc=loc2;
    CCmax=CC(rloc,cloc); 
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n); 
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax)); 
    md2 = fix(m/2); 
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end

    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    output=[error,diffphase,row_shift,col_shift];
    
% Partial-pixel shift
else
    
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data in a 2x larger array
    [m,n]=size(buf1ft);
    mlarge=m*2;
    nlarge=n*2;
    CC=zeros(mlarge,nlarge);
    CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));
  
    % Compute crosscorrelation and locate the peak 
    CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    CCmax=CC(rloc,cloc);
    
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak 
    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2 
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    row_shift=row_shift/2;
    col_shift=col_shift/2;

    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2,
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac; 
        col_shift = round(col_shift*usfac)/usfac;     
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
        % Locate maximum and map back to original pixel grid 
        [max1,loc1] = max(CC);   
        [max2,loc2] = max(max1); 
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);
        rg00 = dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
        rf00 = dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);  
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;    

    % If upsampling = 2, no additional pixel shift refinement
    else    
        rg00 = sum(sum( buf1ft.*conj(buf1ft) ))/m/n;
        rf00 = sum(sum( buf2ft.*conj(buf2ft) ))/m/n;
    end
    error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    % If its only one row or column the shift along that dimension has no
    % effect. We set to zero.
    if md2 == 1,
        row_shift = 0;
    end
    if nd2 == 1,
        col_shift = 0;
    end
    output=[error,diffphase,row_shift,col_shift];
end  

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0),
    [nr,nc]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(i*diffphase);
end
return


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over play.
function play_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%play_Callback(hObject, eventdata, handles)
disp('hey')


% --- Executes on key press with focus on play and none of its controls.
function play_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


function updateimage(cind, handles)
%    guidata(hObject, handles);
    global red 
    global green
    v=get(handles.popupmenu1,'Value'); %% Channel
    x=handles.nx;
    y=handles.ny;
    handles;
    usemean = get(handles.meanimage,'Value');
    handles.avgOn=1;
     if v==1,
        if usemean ==1
            raw = double(handles.gmImg2);
        elseif usemean == 2
            raw=double(green{cind});
          elseif usemean == 3
            cinds = (-3:4) + cind;
            % Boundary conditions
            if cinds(1) <1
                cinds = (cinds - cinds(1)) + 1;  
            end
            if cinds(end)>handles.length
                cinds = cinds - (cinds(end) - handles.length);
            end
            raw = double(green{cinds(1)});
            for j = cinds(2:end)
                raw = raw + double(green{j});
            end
            raw = raw./8;
        end
      
        if (get(handles.checkbox1,'Value'))
          gmin = min(raw(:));
          gmax = max(raw(:));
          set(handles.gmin,'String',num2str(round(gmin)));
          set(handles.gmax,'String',num2str(round(gmax)));
        else
          gmin = (str2num(get(handles.gmin,'String')));
          gmax = (str2num(get(handles.gmax,'String')));
        end
%        disp(round([min(raw(:)) max(raw(:))]))
        
        raw=raw-gmin;
        raw=raw/(gmax-gmin);
%        disp(round([min(raw(:)) max(raw(:))]))
        raw(find(raw<0))=0;
        raw(find(raw>1))=1;
 
        img=cat(3,raw,raw,raw);
    elseif v==2,
        if usemean == 1
            raw = double(handles.rmImg2);
        elseif usemean == 2
            raw=double(red{cind});
        elseif usemean == 3
            cinds = (-3:4) + cind;
            % Boundary conditions
            if cinds(1) <1
                cinds = (cinds - cinds(1)) + 1;  
            end
            if cinds(end)>handles.length
                cinds = cinds - (cinds(end) - handles.length);
            end
            raw = double(red{cinds(1)});
            for j = cinds(2:end)
                raw = raw + double(red{j});
            end
            raw = raw./8;
        end
        
        if (get(handles.checkbox1,'Value'))
          rmin = min(raw(:));
          rmax = max(raw(:));
          set(handles.rmin,'String',num2str(round(rmin)));
          set(handles.rmax,'String',num2str(round(rmax)));
        else
          rmin = (str2num(get(handles.rmin,'String')));
          rmax = (str2num(get(handles.rmax,'String')));
        end
        
        raw=raw-rmin;
        raw=raw/(rmax-rmin);
        raw(find(raw<0))=0;
        raw(find(raw>1))=1;
        img=cat(3,raw,raw,raw);
    else
        if usemean ==1
            raw1 = double(handles.gmImg2);
            raw2 = double(handles.rmImg2);
        elseif (usemean==0) %(usemean==2) || (usemean==3)
            raw1=double(green{cind});
            raw2=double(red{cind});
        end
        if (get(handles.checkbox1,'Value'))
          gmin = min(raw1(:));
          gmax = max(raw1(:));
          set(handles.gmin,'String',num2str(round(gmin)));
          set(handles.gmax,'String',num2str(round(gmax)));
        else
          gmin = (str2num(get(handles.gmin,'String')));
          gmax = (str2num(get(handles.gmax,'String')));
        end
        
        raw1=raw1-gmin;
        raw1=raw1/(gmax-gmin);
        raw1(find(raw1<0))=0;
        raw1(find(raw1>1))=1;

        if (get(handles.checkbox1,'Value'))
          rmin = min(raw2(:));
          rmax = max(raw2(:));
          set(handles.rmin,'String',num2str(round(rmin)));
          set(handles.rmax,'String',num2str(round(rmax)));
        else
          rmin = (str2num(get(handles.rmin,'String')));
          rmax = (str2num(get(handles.rmax,'String')));
        end
        raw2=raw2-rmin;
        raw2=raw2/(rmax-rmin);
        raw2(find(raw2<0))=0;
        raw2(find(raw2>1))=1;
        
        img=cat(3,raw2,raw1,zeros(x,y));
    end       
    set(handles.rimg,'CData',img);
    drawnow;
    
    if 0
    if (usemean == 3)  % Reset useman dialog
      set(handles.meanimage,'Value',2);
    end
    end
    % update stim number
    
    stimwin = floor(4 * (7910/256));
    
    if isfield(handles,'stim') && length(handles.stim>1)
      pind = find(handles.stim(:,1)<cind,1,'last');
      if ~isempty(pind) && ((cind-handles.stim(pind,1))<stimwin)
        set(handles.stimtext,'String',num2str(handles.stim(pind,2)),'ForegroundColor',[1 0 0]);
      else
        set(handles.stimtext,'String','','ForegroundColor',[0 0 0]);
      end
    else
        set(handles.stimtext,'String','','ForegroundColor',[0 0 0]);
    end
    

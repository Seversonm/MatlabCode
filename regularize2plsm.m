function success = regularize2plsm(fname,chan)
% function success = regularize2plsm(fname,chan)
% chan = 1 red
% chan = 2 green
% chan = 0 or anything else, no regularization

  if ~exist('chan')
      chan = 1;
  end

  zout = fliplr([45 47 124 92 ]);

  framefinder(fname,0);
  fid = fopen([fname '_FF'],'r');
  a = fread(fid,'uint64');
  framelist= reshape(a,2,size(a,1)/2)';
  fclose(fid);
  
  
  
  fid = fopen(fname,'r','b');
  if fid<0
   success = 0;
   return;
  end

  [desc, gccode,fid2,red,green,vdat] = l2phfile_resFFJ(fid,0,1,framelist);
  mgimg=green{1};
  mrimg = red{1};
  green{1} = int16(green{1});
  red{1} = int16(red{1});
  offsets = zeros(size(framelist,1),2);
  fprintf('\n');
  for i=2:size(framelist,1)
%i
    [desc, gcc,fid2,r1,g1,vd] = l2phfile_resFFJ(fid,0,i,framelist);
    if chan==2
      [p]=dftregistration(fft2(mgimg),fft2(double(g1{1})),10);
      green{i}=int16(real(ifft2(dftshift(fft2(double(g1{1})),p(2),p(3),p(end)))));
      red{i}=int16(real(ifft2(dftshift(fft2(double(r1{1})),p(2),p(3),p(end)))));
    elseif chan==1
      [p]=dftregistration(fft2(mrimg),fft2(double(r1{1})),10);
      green{i}=int16(real(ifft2(dftshift(fft2(double(g1{1})),p(2),p(3),p(end)))));
      red{i}=int16(real(ifft2(dftshift(fft2(double(r1{1})),p(2),p(3),p(end)))));
    else
      green{i} = g1{1};
      red{i} = r1{1};
    end
%    green{i}=int16(real(ifft2(dftshift(fft2(double(g1{1})),p(2),p(3),p(end)))));
%    red{i}=int16(real(ifft2(dftshift(fft2(double(r1{1})),p(2),p(3),p(end)))));
    gccode{i} = gcc{1};
    vdat{i} = vd{1};
    mgimg = mgimg+double(green{i});
    mrimg = mrimg + double(red{i});
    offsets(i,:) = offsets(i,:) + p(3:4);
   fprintf('%c',[8 zout(rem(i,4)+1)])
 
  end
  fclose(fid)

 fnout = [fname(1:end-5) 'mat'];
 save(fnout,'desc','gccode','green','red','vdat','offsets','-v7.3');


function [] = framefinder(fname,verbose)
  
if ~exist('verbose')
  verbose = 0;
end

% If it already exists, return
fid = fopen([fname '_FF'],'r','b');
if fid>2
  return;
end


fid = fopen([fname],'r','b');

if fid<0
  desc =[]; rtimes = []; dualphdata  = [];
  return;
end
fseek(fid,0,1);
endoffile = ftell(fid);
fseek(fid,0,-1);

cnt = 1;
done = 0;
fdat = [];
while ~done
  fdat =[fdat; cnt ftell(fid)];

  bnel = fread(fid,1,'uint32');
  for j=1:bnel
    nel = fread(fid,1,'uint32');
    fseek(fid,nel,'cof');
  end
  if ftell(fid)==endoffile
    done = 1;
  end
  cnt = cnt + 1;
end
fclose(fid)

fid = fopen([fname '_FF'],'w')
fwrite(fid,fdat','uint64')
%fprintf(fid,'%12g %12g\n',fdat')

fclose(fid)




function [desc,gccode,fid,red,green,vdat] = l2phfile_resFFJ(fid,verbose,frnum, framelist)
% function [desc,rtimes,fid,vdat] = l2phfile(pathname,fname,verbose,frnum)
if ~exist('verbose')
  verbose = 0;
end

framelocs = framelist(frnum,2);
fseek(fid,0,1); 
endoffile = ftell(fid);

framecnt=0;
templatecnt=0;

cnt = 1;
done = 0;
first = 1;
lframe= 1;

for cframe = 1:length(frnum) 
    
%  disp('hey')   
  fseek(fid,framelocs(cframe),-1);

  bnel = fread(fid,1,'uint32');
% Element 1
  dat = unflatten2string(fid,'char',verbose,lframe);
% Element 2
  tim = unflatten2string(fid,'char',verbose,lframe);
% Element 3
  junk = unflatten2string(fid,'char',verbose,lframe);

% Element 4
  nbytes = fread(fid,1,'uint32');
  nel = fread(fid,1,'uint32');

  parnames = [];
  for j=1:nel
    parnames{j} = unflatten2string(fid,'char',verbose,lframe);
  end

% Element 5
  parvals = unflatten2string(fid,'1darraydouble',verbose,lframe);


  if lframe
    if first
      d1 = [];
      for j=1:length(parvals)
        d =setfield(d1,rmwhitespace(parnames{j}),parvals(j));
        d1 = d;
      end
      desc(cframe) = d1;
    else
      desc(cnt).AOrate = 0;
      for j=1:length(parvals)
        desc(cframe)=setfield(desc(cnt),rmwhitespace(parnames{j}),parvals(j));
      end
    end
    if first
      oldvals =[desc(cframe).Xres desc(cframe).Yres desc(cframe).Linescan];
    end
    if sum(oldvals - [desc(cframe).Xres desc(cframe).Yres desc(cframe).Linescan]) > 0
      disp(['CHANGED DATA TYPE at header ' num2str(cframe) ':  ' num2str(oldvals) ' -> ' num2str([desc(cframe).Xres desc(cframe).Yres desc(cframe).Linescan])])
    end
    first =0;
  end
  
  
 % Element 6
 rr = unflatten2string(fid,'char',verbose,lframe);
 
 % Element 7
 vdat{cframe} = unflatten2string(fid,'2darray',verbose,lframe);
 
 % Element 8
 comment{cframe} = unflatten2string(fid,'char',verbose,lframe);
 % Element 9
 %dont want in this data for bulk loading
 
%  red{cframe} = double(unflatten2string(fid,'2darray',verbose,lframe));

dat = unflatten2string(fid,'2darray',verbose,lframe);
green{cframe} = double(dat);
 
 % Element 10
 dat=unflatten2string(fid,'2darray',verbose,lframe);
 red{cframe} = double(dat);%double(dat); %%%%%%%%%%%%%%%%%%%%%%

  if bnel>10
      % Element 11
    gccode{cframe} = unflatten2string(fid,'char',verbose,lframe);
    %disp(gccode{cframe}); fprintf('\n');
    for cel = 12:bnel
      unflatten2string(fid,'unknown',verbose,lframe);
    end
  end
end % while ~done

desc = desc(1);


function strout = rmwhitespace(strin)
  strout= strin;
  whitespace = [9 10 11 12 13 32];
  for j=1:length(whitespace)
    z1 = find(strout==char(whitespace(j)));
    strout(z1) = [];
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
    disp('H')
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(i*diffphase);
    disp('I')
end
return



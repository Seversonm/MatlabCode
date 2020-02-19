function [dat] = unflatten2string(fid,strtype,verbose,doload)
% function [dat] = unflatten2string(fid,strtype,verbose,doload)
dat = [];
nel = fread(fid,1,'uint32');

if ~doload
  fseek(fid,nel,'cof');
  return;
end

if verbose
  disp(['Reading ' num2str(nel) ' bytes'])
end

if strmatch(strtype,'char')
  dat = char(fread(fid,nel,'char'))';
  return;
end

if strmatch(strtype,'unknown')
  fseek(fid,nel,'cof');
end

if strmatch(strtype,'double')
  dat = (fread(fid,nel,'double'))';
  return;
end

%

if strmatch(strtype,'float')
  dat = (fread(fid,nel,'float32'))';
  return;
end


if strmatch(strtype,'uint16')
  dat = (fread(fid,nel,'uint16'))';
  return;
end

if strmatch(strtype,'int16')
  dat = (fread(fid,nel,'int16'))';
  return;
end

if strmatch(strtype,'uint32')
  dat = (fread(fid,nel,'uint32'))';
  return;
end

if strmatch(strtype,'int32')
  dat = (fread(fid,nel,'int32'))';
  return;
end

if strmatch(strtype,'2darray')
  xlen = fread(fid,1,'uint32');
  ylen = fread(fid,1,'uint32');
  dat2 = (fread(fid,xlen*ylen,'int16'))';
  dat = zeros(xlen,ylen,'int16');
  for j=1:xlen
    dat(j,1:ylen) = dat2((1:ylen) + (j-1)*ylen);
  end
  return;
end

if strmatch(strtype,'1darray')
  xlen = fread(fid,1,'uint32');
  dat = (fread(fid,xlen,'int16'))';
  return;
end


if strmatch(strtype,'1darrayfloat')
  xlen = fread(fid,1,'int32');
  dat = (fread(fid,xlen,'float32'))';
  return;
end

if strmatch(strtype,'1darraydouble')
  xlen = fread(fid,1,'int32');
  dat = (fread(fid,xlen,'float64'))';
  return;
end

end

  

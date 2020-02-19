function retvec = resamplevec(vec,samplerate,int)
% function retvec = resamplevec(vec,samplerate,int)
%    int = 1, try to linear interpolate the vector
%    int = 0, no interpolation


if ~exist('int')
  int = 1;
end
if (samplerate > 1)
if int
  retvec = zeros((length(vec)-1) * samplerate + 1 ,1);
  for j=1:length(vec) - 1
    step = (vec(j+1) - vec(j))/samplerate;
    for  k=1:samplerate
      retvec((j-1)*samplerate + k ) = vec(j) + (k-1) * step;
    end
  end;
  retvec((j) * samplerate + 1) = vec(length(vec));
else
  retvec = zeros(length(vec)*samplerate,1); 
  for j=1:length(vec)
    retvec((j-1)*samplerate+1:j*samplerate) = ones(samplerate,1)*vec(j);
  end
end
else   % if samplerate is less than one .. need to take an average
       % This is built to only rebin according to 0.5, 0.25 0.125, ...

  binsize = 1/samplerate;
  retvec = zeros(ceil(length(vec) * samplerate),1);
  cnt = 1;
  for j=1:binsize:length(vec)
    retvec(cnt) = mean(vec(floor(j):min([floor(j+binsize - 1) length(vec)])));
    cnt = cnt + 1;
  end
end

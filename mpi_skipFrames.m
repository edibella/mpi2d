function cinemri1=mpi_skipFrames(framesToSkip, dynamicData, ranget)
% framesToSkip is vector of frame numbers to skip

if nargin==0,
    error('need arguments in mpi_skipFrames');
end

%   if framesToSkip(1)==1
%       framesToSkip=framesToSkip(2:end);
%   end
  
%   ii=find(framesToSkip>size(dynamicData,3))
%   framesToSkip=framesToSkip(1:min(ii)-1);
%   % need ranget
  % keep only frames within ranget:
  ii=find(framesToSkip< min(ranget))
  if ~isempty(ii)
      framesToSkip=framesToSkip(max(ii)+1:end);
  end
  ii=find(framesToSkip> max(ranget)-1)
  if ~isempty(ii)
      framesToSkip=framesToSkip(1:min(ii)-1);
  end
  framesToSkip= framesToSkip - min(ranget) +1;
  
  cinemri1=dynamicData;

  nSetsOfFrames=1;
  nFramesInSet=zeros(1,512);   % was 100  8/11  EVRD

  for n=1:length(framesToSkip)-1
     if (framesToSkip(n) == framesToSkip(n+1)-1)
        nFramesInSet(nSetsOfFrames)=nFramesInSet(nSetsOfFrames)+1;
        if n==length(framesToSkip)-1
           nFramesInSet(nSetsOfFrames)=nFramesInSet(nSetsOfFrames)+1;
        end
     else
        nFramesInSet(nSetsOfFrames)=nFramesInSet(nSetsOfFrames)+1;
        nSetsOfFrames=nSetsOfFrames+1;
        if n==length(framesToSkip)-1
           nFramesInSet(nSetsOfFrames)=nFramesInSet(nSetsOfFrames)+1;
        end
     end
  end

  %handle first set differently if starts at 1. Will replace it just with
  %rhs.. 
  
  
% for each set, do linear interp if nFramesInSet is 1
counter=1;
for c=1:nSetsOfFrames
  if framesToSkip(counter)==1
     rhsValue=squeeze(dynamicData(:,:,framesToSkip(counter)+nFramesInSet(c)));
     for m=1:nFramesInSet(c)
        cinemri1(:,:,framesToSkip(counter))= rhsValue;
        counter=counter+1;
     end
  else
      
  lhsValue=squeeze(dynamicData(:,:,framesToSkip(counter)-1));
  rhsValue=squeeze(dynamicData(:,:,framesToSkip(counter)+nFramesInSet(c)));
  for m=1:nFramesInSet(c)
     cinemri1(:,:,framesToSkip(counter))= lhsValue + (m/(nFramesInSet(c)+1))*(rhsValue - lhsValue);
     cinemri1(:,:,framesToSkip(counter))==0;
     counter=counter+1;
  end
  end
  
end

return

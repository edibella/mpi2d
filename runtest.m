function  numRuns=runtest(err)
% return number of positive or neg. runs in vector err

nTimes=length(err);

numRuns=0;
for ii=1:nTimes-1
   if sign(err(ii)) ~= sign(err(ii+1))    % next value same sign as last value 
      numRuns=numRuns+1;
   end
end

return

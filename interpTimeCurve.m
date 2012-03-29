function  newTimeCurve=interpTimeCurve(selectedTimeStamps, delta_t, origTimeCurve)
% newTimeCurve will have length totaltime/delta_t
% example call: 
%   bldcurve=interpTimeCurves(timeStampFileAIF, delta_t, bldcurve);

showfigs=0;


%    load(timeStampFile);   % assumes gives variable timeStamp
%    if(size(timeStamp,1) == 1)
%        %the timestamp was not saved in the format mpi2d expects
%        timeStamp = timeStamp';
%        save(timeStampFile,'timeStamp');
%    end
   timeStamps=squeeze(selectedTimeStamps(:,1)-selectedTimeStamps(1,1));   % slice 1 assumedly
   %keyboard
   totalTime=selectedTimeStamps(length(origTimeCurve),1)-selectedTimeStamps(1,1);
   nTimes=round(totalTime/delta_t);

   % units are in seconds after subtract off first time
   origSampleTimes= timeStamps;
   if(length(origSampleTimes) > length(origTimeCurve))
       origSampleTimes=origSampleTimes(2:length(origTimeCurve)+1);
   end

if (showfigs)
   figure; clf; hold on
   set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
   plot((1:length(origTimeCurve))*max(origSampleTimes)/(length(origTimeCurve)),origTimeCurve,'-or')
   plot(origSampleTimes,origTimeCurve,'-xm')
   legend('assume uniform','with time stamps')
end

%   extrapval=mean(bldcurve(nTimes-4:nTimes-1));
   extrapval=0;

   newTimeCurve = interp1(origSampleTimes, origTimeCurve, delta_t*(0:nTimes-1),'cubic',extrapval);
% NOTE: interp1 seems different in matlab6! Using matlab7 here!
   newTimeCurve=newTimeCurve';
%   for zz=1:10
%        if bldcurve(zz)==extrapval
%            bldcurve(zz)=0;
%        end
%   end

if (showfigs)   % this one is if have dual bolus.  Need to hack if want
%  to see it... Maybe copy into fit_*.m
   set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
   bb=delta_t*(0:nTimes-1); bb=bb';
%   cc=bb-origSampleTimes0;
 %  dd=bb-origSampleTimes;
%   figure(9)
 %  clf; hold on;
%   plot(cc+origSampleTimes0(1),'g','linewidth',3)
%   plot(dd+origSampleTimes(1),'b','linewidth',3)
%   legend('1st bolus', '2nd bolus')
%   disp('delta_t is'), delta_t
%   ylabel('Seconds')
%   xlabel('Frame number')

%    figure(1); clf; hold on
%   plot(delta_t*(0:nTimes-1),origbldcurve(:),'-dg')
%   plot(delta_t*(0:new_nTimes-1),bldcurve(:),'-dk')
%    plot(delta_t*(0:nTimes-1),origtisscurve(:,3),'-dg')
%    plot(delta_t*(0:new_nTimes-1),tisscurve(:,3),'-dm')
%    plot(origSampleTimes,origbldcurve,'-xr')
%plot(origSampleTimes,origtisscurve(:,3),'-xk')
%    xlabel('seconds')

    figure(20); clf;
    set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
    plot(origSampleTimes(2:end) - origSampleTimes(1:end-1),'-om','linewidth',2.2)
    xlabel('Frame number'); ylabel('Interframe interval (seconds)')
end

%disp('showing for impact of time stamps')

return
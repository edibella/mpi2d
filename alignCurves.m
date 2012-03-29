function  [t_delay, curveShifted]=alignCurves(bldcurve, sat_bldcurve)
% shifts first curve to match second in time. Only by integers, uses upslopes (or cross-corr. if chg)

numSkipUpslope=3;

nTimes=length(bldcurve);

% clf; figure(1);
% plot(bldcurve)
% hold on
% plot(sat_bldcurve,'r')

% 
[aa]=xcorr(bldcurve, sat_bldcurve);
[maxval indexShift]=max(aa);
timeDiff_AIF=nTimes-indexShift;


%call upslope routine:

N=3;    % Schwitter2001 used 3 point linear fit in blood and 5 point in tissue
N=4; % changed 10/23/05 to be more general? we're using every beat...
maxslope_bld=0; xvec=1:N; xvec=xvec';
%slope_location_bld=1;
for i=numSkipUpslope:nTimes-N
   bb = find(bldcurve(i:i+N-1,1)~=0);
   if (~isempty(bb))
      [tmpcoeffs,Serr] = lregress(bldcurve(i:i+N-1,1),xvec);
      slope=tmpcoeffs(2,1);
   else
      slope=0;
   end
   if slope > maxslope_bld
       maxslope_bld = slope;
       slope_location_bld=i;   % can also use Serr to check if decent fit
   end
end
N=3;    % Schwitter2001 used 3 point linear fit in blood and 5 point in tissue
N=6; % changed 10/23/05 to be more general? we're using every beat...
N=4; % should match if aligning bld curves!!! Can be longer if aligning tiss and bld
maxslope_bld=0; xvec=1:N; xvec=xvec';
%slope_location_sat=1;
for i=numSkipUpslope:nTimes-N
   bb = find(sat_bldcurve(i:i+N-1,1)~=0);
   if (~isempty(bb))
      [tmpcoeffs,Serr] = lregress(sat_bldcurve(i:i+N-1,1),xvec);
      slope=tmpcoeffs(2,1);
   else
      slope=0;
   end
   if slope > maxslope_bld
       maxslope_bld = slope;
       slope_location_sat=i;   % can also use Serr to check if decent fit
   end
end

% clf; hold on
%    plot(1:nTimes,bldcurve(:,1),'r','Linewidth',2.5)
%    tmpcoeffs = mpi_twoplot(bldcurve(slope_location_bld:slope_location_bld+N-1),xvec,slope_location_bld); % to see plots
%    plot(1:nTimes,sat_bldcurve(:,1),'k','Linewidth',2.5)
%    tmpcoeffs = mpi_twoplot(sat_bldcurve(slope_location_sat:slope_location_sat+N-1),xvec,slope_location_sat); % to see plots
% 

t_delay=-timeDiff_AIF;
disp('above delay is from cross correlation, using delay below from upslopes')
t_delay=slope_location_bld -slope_location_sat;

tmpbldcurve=bldcurve;
if t_delay>0
  disp('problem with shifting non-sat. AIF to left. Havent implemented that yet') 
% hmm, should be ok for integer shifts! Test again.
    tmpbldcurve(1:nTimes-fix(t_delay))=bldcurve(1+fix(t_delay):nTimes);
    tmpbldcurve(nTimes:-1:nTimes-fix(t_delay))=0;
    t_delay=t_delay-fix(t_delay);
    bldcurveShifted = tmpbldcurve;   %interp1(1:nTimes,tmpbldcurve,1+t_delay:nTimes+t_delay,'splines',0);
else
    tmpbldcurve(1-fix(t_delay):nTimes)=bldcurve(1:nTimes+fix(t_delay));
    tmpbldcurve(1:1-fix(t_delay))=0;
    t_delay=t_delay-fix(t_delay);
    bldcurveShifted = interp1(1:nTimes,tmpbldcurve,1+t_delay:nTimes+t_delay,'splines',0);
end

figure(2); clf; hold on
   plot(1:nTimes,bldcurve(:,1),'r','Linewidth',2.5)
   plot(1:nTimes,bldcurveShifted,'g','Linewidth',2.5)
   plot(1:nTimes,sat_bldcurve(:,1),'k','Linewidth',2.5)
   legend('firstcurve','firstcurveShifted','satbldcurve or 2nd curve')
drawnow

t_delay=slope_location_bld -slope_location_sat;
curveShifted=bldcurveShifted;
return;

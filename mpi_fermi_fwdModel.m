%  mainly so can change how create fittedcurve  in only one place!!!

function [fittedcurve, tmpbldcurve]=mpi_fermi_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi, kwo, t_delay, T0, est_fb, spillover);   % note all parameters are scalars!!!

nTimes=length(bldcurve);
tmp=0:nTimes-1;
ttimes = (delta_t/60)*tmp';
imp=zeros(nTimes,1);
imp(1)=1.0;

F=kwi;

if ~exist('est_fb')
  est_fb=0;
end
if ~exist('spillover')
  spillover=0;
end


% adding 7/2/04 to look at - even small say .05 fv cuts k1 in half..
%est_fb=0.07;
%t_delay=t_delay*0;
t_delayVp=t_delay;

% t_delay just in the blood curve currently 4/04
% is t_delay always negative?  interp1 just doesn't give points < 0
    tmpbldcurve=bldcurve;
    tmpbldcurve(1-fix(t_delay):nTimes)=bldcurve(1:nTimes+fix(t_delay));
    tmpbldcurve(1:1-fix(t_delay))=0;
    t_delay=t_delay-fix(t_delay);
%    tmpbldcurve = interp1(1:nTimes,tmpbldcurve,1+t_delay:nTimes+t_delay,'splines',0);
    tmpbldcurve = interp1(1:nTimes,tmpbldcurve,1+t_delay:nTimes+t_delay,'cubic',0);
% changed 9/19/05 to consider est_fb just models spillover. So no delay in it!
%    fittedcurve = conv(tmpbldcurve,(delta_t/60*F./(exp(delta_t/60*kwo*(ttimes-T0))+1) ));
% convert T0 to seconds rather than frames, 10/16/05
%    fittedcurve = conv(tmpbldcurve,((1.0/60)*F./(exp(kwo*(ttimes-T0))+1) ));

% 3/27/06  changed back to delta_t in front
%   fittedcurve = conv(tmpbldcurve,((delta_t/60)*F./(exp(kwo*(ttimes-T0))+1) ));
% 5/19/06, make T0 in seconds
    fittedcurve = conv(tmpbldcurve,((delta_t/60)*F./(exp(kwo*(ttimes-T0/60.0))+1) ));
    
    h=((delta_t/60)*F./(exp(kwo*(ttimes-T0/60.0))+1) ); save hTMP.mat h est_fb
        
    tmpsat_bldcurve=sat_bldcurve;
    tmpsat_bldcurve(1-fix(t_delayVp):nTimes)=sat_bldcurve(1:nTimes+fix(t_delayVp));
    tmpsat_bldcurve(1:1-fix(t_delayVp))=0;
    t_delayVp=t_delayVp-fix(t_delayVp);
%    tmpsat_bldcurve = interp1(1:nTimes,tmpsat_bldcurve,1+t_delayVp:nTimes+t_delayVp,'splines',0);
    tmpsat_bldcurve = interp1(1:nTimes,tmpsat_bldcurve,1+t_delayVp:nTimes+t_delayVp,'cubic',0);
    fittedcurve=fittedcurve(1:nTimes);
    
    fittedcurve=(1-spillover)*fittedcurve(1:nTimes) + est_fb*tmpsat_bldcurve(1:nTimes)' + spillover*sat_bldcurve(1:nTimes);
    
% trying "bad idea" again 3/29/06
    

% made match with sat*fwd*.m  4/4/06


%figure(2); clf;
%plot(fittedcurve)
%figure(3); clf;
%impulseResponse=delta_t/60*F./(exp(delta_t/60*kwo*(ttimes-T0))+1);
%impulseResponse(T0:end)=0
%plot(impulseResponse)
%plot(delta_t/60*F./(exp(delta_t/60*kwo*(ttimes-T0))+1))
%%hold on
%plot(delta_t/60*F./(exp(delta_t/60*kwo*(ttimes))+1),'r')
%drawnow

return

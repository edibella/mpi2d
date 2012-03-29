%  mainly so can change how create fittedcurve  in only one place!!!

function [fittedcurve, tmpbldcurve]=mpi_PET_fwdModel(bldcurve, sat_bldcurve, delta_t, ptf, mbf, t_delay, est_fb, spillover);   % note all parameters are scalars!!!
%function [fittedcurve, tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi, PartCoeff, t_delay, est_fb, spillover);   % note all parameters are scalars!!!


% C(t)=ptf*kwi*conv(bld,exp(-kwi/.91)) + vp*bld  + vp_rv*bldRV

nTimes=length(bldcurve);
tmp=0:nTimes-1;
ttimes = tmp';
imp=zeros(nTimes,1);
imp(1)=1.0;

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
%    fittedcurve = conv(tmpbldcurve,((1.0/60)*kwi*exp(-kwo*delta_t/60*(ttimes)) ) );
% changed 3/15/06 to normalize for more samples in conv
    
    fittedcurve = conv(tmpbldcurve,((delta_t/60)*ptf*mbf*exp(-(mbf./0.91)*delta_t/60*(ttimes)) ) );
%    fittedcurve = conv(tmpbldcurve,((delta_t/60)*kwi*exp(-(kwi/PartCoeff)*delta_t/60*(ttimes)) ) ); % changed to directly estimate the Partition Coefficient (NP 10/10/06)
  
%    fittedcurve = conv(tmpbldcurve,((1-spillover)/60*kwi*exp(-kwo*delta_t/60*(ttimes)) ) );  % this was a bad idea. RIP 10/18/05
    tmpsat_bldcurve=sat_bldcurve;
    tmpsat_bldcurve(1-fix(t_delayVp):nTimes)=sat_bldcurve(1:nTimes+fix(t_delayVp));
    tmpsat_bldcurve(1:1-fix(t_delayVp))=0;
    t_delayVp=t_delayVp-fix(t_delayVp);
    tmpsat_bldcurve = interp1(1:nTimes,tmpsat_bldcurve,1+t_delayVp:nTimes+t_delayVp,'cubic',0);
    fittedcurve=fittedcurve(1:nTimes);
    tmpfittedcurve=fittedcurve(1:nTimes);
    spillover=0;
    
   % fittedcurve=(1-spillover)*fittedcurve(1:nTimes) + est_fb*tmpsat_bldcurve(1:nTimes)' + spillover*sat_bldcurve(1:nTimes);
% trying "bad idea" again 3/29/06


fittedcurve=(1-spillover)*fittedcurve(1:nTimes) + est_fb*tmpbldcurve(1:nTimes)' + spillover*sat_bldcurve(1:nTimes); %EVRD 8/27/11, seems
% shoould use same AIF for vb and AIF..  If saturated in blood, that is
% like assuming water exchange slow between capillary and tissue for
% capillaries.  Hmm, which Buckley said should look at more. But most work
% looks at in and out of cells, within the interstitium. 

% clf; hold on
%      set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
% %      plot(ttimes,bldcurve,'r','LineWidth',3)
%      plot(ttimes,fittedcurve,'-','LineWidth',2.2)
%      plot(ttimes,tmpfittedcurve,'x','LineWidth',2.2)
%      plot(ttimes,est_fb*tmpsat_bldcurve(1:nTimes)','m','LineWidth',2)
%      plot(ttimes,spillover*sat_bldcurve(1:nTimes),'c','LineWidth',2)
% 
%      legend('Total Fit','non-blood part', 'Vb part', 'Spillover part')
%       xlabel('Time (sec)')
% %      ylabel('Gd concentration or deltaSI')
%      ylabel('deltaSI')
% drawnow
% t_delay
% pause


return;

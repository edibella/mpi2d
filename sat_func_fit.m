function [Fsum]=sat_func_fit(x,bldcurve, sat_bldcurve, tisscurve,delta_t,nRegs,fixedDelay, fixedVp, estGlobalDelay)
%  matlab function to be called by fmincon


nTimes=length(bldcurve);

tmp=0:nTimes-1;
ttimes = tmp';

kwi = x(1:nRegs);
kwo = x(nRegs+1:2*nRegs);
%PartCoeff = x(nRegs+1:2*nRegs); % this was added to directly estimate the Partition Coefficient instead of kwo above (NP 10/10/06)

t_delay = fixedDelay*ones(1,nRegs);
est_fb = fixedVp*ones(1,nRegs);
spillover = zeros(1,nRegs);
if (length(x)>2*nRegs)
   t_delay = fixedDelay*ones(1,nRegs);
   if fixedDelay==99
       t_delay = x(2*nRegs+1:3*nRegs);
   else
       est_fb = x(2*nRegs+1:3*nRegs);
   end
   if fixedDelay~=99 && fixedVp~=99
       spillover = x(2*nRegs+1:3*nRegs);
   end
end

if (length(x)>3*nRegs)
   if fixedDelay==99 && fixedVp==99
       est_fb = x(3*nRegs+1:4*nRegs);
   else
       spillover = x(3*nRegs+1:4*nRegs);
   end
end

if (length(x)>4*nRegs)  % means delay, est_fb, and spillover being estimated
       spillover = x(4*nRegs+1:5*nRegs);
end

% test all combinatiosn of fixedDelay, fixedVp, and size of x
%keyboard
%%%%%%%% CHANGED NP 071807 to force deltaT and Vb to be zero (also see fit_model_2comp.m) 
% est_fb=0;
% t_delayVp=0;
% spillover=0;

%t_delayVp = x(4*nRegs+1:5*nRegs);
% t_delayVp=0;
% spillover=0;


jitter=0; newRate=1;    % for downsampling ISMRM 2003 (not used here)

imp=zeros(nTimes,1);
imp(1)=1.0;                                  

errsum = 0.0;

%load Output/deltaSIcurves.study8.slice3.AIF_8_3_1.mat
%bldcurve=deltaSIcurves(1,:);
%sat_bldcurve=bldcurve';
for i=1:nRegs
  [est_curves, tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), spillover(i));
%  [est_curves, tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), PartCoeff(i), t_delay(i), est_fb(i), spillover(i)); % changed to directly estimate the Partition Coefficient (NP 10/10/06)

%  [est_curves, tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), spillover(i));
%  [est_curves, tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), t_delayVp(i), spillover(i));
%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'linear',0);
% try this 4/1/04
%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'v5cubic',0);
%%%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'splines',0);
%      est_curves = conv(tmpbldcurve, delta_t/60*((1-est_fb(i))*kwi(i)*exp(-kwo(i)*delta_t/60*(ttimes-t_delay(i)*ones(nTimes,1))) + est_fb(i)*imp) );
%%%      est_curves = conv(tmpbldcurve, ((1-est_fb(i))*(delta_t/60)*kwi(i)*exp(-kwo(i)*delta_t/60*(ttimes)) + est_fb(i)*imp) );
      
est_curves=est_curves((jitter+1):newRate:nTimes);   % downsampling step
     err = (tisscurve(1:nTimes,i) - est_curves);

%     load Output/deltaSIcurves.study998.slice3.AIF_998_3_1.mat
%     load Output/deltaSIcurves.study8.slice3.AIF_8_3_1.mat
%aa=bldcurve - deltaSIcurves(1,:)
%keyboard
      
     errsum = errsum + sum(err.*err);

if estGlobalDelay~=99 && fixedDelay==99
     lambda=.33*errsum;  % changed from .1 11/1/05, from .2 11/11
                        % changed from .5 to .33 4/22/06
     lambda=0;   % changed 5/3/06
     errsum = errsum + lambda*(t_delay(i)-estGlobalDelay)^2;
end

% clf 
% hold on
% plot(bldcurve,'r')
% plot(tmpbldcurve,'m')
% plot(est_curves(1:nTimes))
% %plot(tisscurve(:,i),'g')
% plot(tisscurve(:,i),'kx')
% drawnow
% pause

end

Fsum = errsum;
%disp('ERROR is '),Fsum
% %t_delay
% clf; plot(est_curves,'k')
% hold on
% plot(tisscurve(:,1),'x')
% drawnow


% constraint:   G < zeros(G)
% constrant a function of x (tissue values) to have norm of one? 
% try just all curves greater than 0 
% or say greater than 1

% not currently using
%G=-10;
% and constrain  fb
%G = -10*est_fb;


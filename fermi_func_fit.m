function [Fsum]=fermi_func_fit(x,bldcurve, sat_bldcurve, tisscurve,delta_t,nRegs,fixedDelay, fixedVp,  estGlobalDelay)
%  matlab function to be called by fmincon

nTimes=length(bldcurve);

% 1/17/06, try shorter curves for Fermi model.
%nTimes=45;
% commented out 3/27/06


tmp=0:nTimes-1;
ttimes = tmp';

kwi = x(1:nRegs);
kwo = x(nRegs+1:2*nRegs);
T0  = x(2*nRegs+1:3*nRegs);

t_delay = fixedDelay*ones(1,nRegs);
est_fb = fixedVp*ones(1,nRegs);
spillover = zeros(1,nRegs);
if (length(x)>3*nRegs)
   t_delay = fixedDelay*ones(1,nRegs);
   if fixedDelay==99
       t_delay = x(3*nRegs+1:4*nRegs);
   else
       est_fb = x(3*nRegs+1:4*nRegs);
   end
   if fixedDelay~=99 && fixedVp~=99
       spillover = x(3*nRegs+1:4*nRegs);
   end
end
if (length(x)>4*nRegs)
   if fixedDelay==99 && fixedVp==99
       est_fb = x(4*nRegs+1:5*nRegs);
   else
       spillover = x(4*nRegs+1:5*nRegs);
   end
end
if (length(x)>5*nRegs)  % means delay, est_fb, and spillover being estimated
       spillover = x(5*nRegs+1:6*nRegs);
end

% hack 5/16/06 to make Vb always zero for Fermi
est_fb = fixedVp*zeros(1,nRegs);

% t_delayVp=0;
spillover=0;

%keyboard

jitter=0; newRate=1;    % for downsampling ISMRM 2003 (not used here)
errsum = 0.0;
for i=1:nRegs
  [est_curves, tmpbldcurve]=mpi_fermi_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), T0(i), est_fb(i), spillover(i));
%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'linear',0);
% try this 4/1/04
%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'v5cubic',0);
%%%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'splines',0);
%      est_curves = conv(tmpbldcurve, delta_t/60*((1-est_fb(i))*kwi(i)*exp(-kwo(i)*delta_t/60*(ttimes-t_delay(i)*ones(nTimes,1))) + est_fb(i)*imp) );
%%%      est_curves = conv(tmpbldcurve, ((1-est_fb(i))*(delta_t/60)*kwi(i)*exp(-kwo(i)*delta_t/60*(ttimes)) + est_fb(i)*imp) );
      
est_curves=est_curves((jitter+1):newRate:nTimes);   % downsampling step
     err = (tisscurve(1:nTimes,i) - est_curves);
      
     %errsum = errsum + sum(err.*err);
     errsum = errsum + sqrt(sum(err.*err)); % using a slightly different erreor metric than above   NP 110308

if estGlobalDelay~=99 && fixedDelay==99
     lambda=.33*errsum;  % changed from .1 11/1/05, from .2 11/11
                        % changed from .5 to .33 4/22/06
     lambda=0;   % changed 5/3/06
     errsum = errsum + lambda*(t_delay(i)-estGlobalDelay)^2;
end


%figure(1); clf; 
%hold on
%plot(bldcurve,'r')
%plot(tmpbldcurve,'m')
%plot(est_curves(1:nTimes))
%plot(tisscurve(:,i),'g')
%plot(tisscurve(:,i),'kx')
%drawnow

%t_delay

end

%Fsum = errsum;
Fsum = errsum/length(err); % here I'm normalizing the errors to the lengths of the measured and estimated data above NP 110308

%disp('ERROR is '),Fsum
%t_delay
%kwi

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


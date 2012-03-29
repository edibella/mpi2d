function [Fsum]=sat_func_fitNoBlood(x,bldcurve, sat_bldcurve, tisscurve,delta_t,nRegs,fixedDelay, fixedVp, estGlobalDelay)
%  matlab function to be called by fmincon


nTimes=length(bldcurve);

tmp=0:nTimes-1;
ttimes = tmp';

kwi = x(1:nRegs);
kwo = x(nRegs+1:2*nRegs);

t_delay = fixedDelay*ones(1,nRegs);
est_fb = fixedVp*ones(1,nRegs);
if (length(x)>2*nRegs)
   t_delay = fixedDelay*ones(1,nRegs);
   if fixedDelay==99
       t_delay = x(2*nRegs+1:3*nRegs);
   else
       est_fb = x(2*nRegs+1:3*nRegs);
   end
end
if (length(x)>3*nRegs)
   est_fb = x(3*nRegs+1:4*nRegs);
end

%t_delayVp = x(4*nRegs+1:5*nRegs);
%spillover = x(5*nRegs+1:6*nRegs);
spillover = x(4*nRegs+1:5*nRegs);

 est_fb=0;
 t_delayVp=0;
 spillover=0;

jitter=0; newRate=1;    % for downsampling ISMRM 2003 (not used here)

imp=zeros(nTimes,1);
imp(1)=1.0;


errsum = 0.0;
for i=1:nRegs
  [est_curves, tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), spillover(i));
%  [est_curves, tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), t_delayVp(i), spillover(i));
%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'linear',0);
% try this 4/1/04
%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'v5cubic',0);
%%%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'splines',0);
%      est_curves = conv(tmpbldcurve, delta_t/60*((1-est_fb(i))*kwi(i)*exp(-kwo(i)*delta_t/60*(ttimes-t_delay(i)*ones(nTimes,1))) + est_fb(i)*imp) );
%%%      est_curves = conv(tmpbldcurve, ((1-est_fb(i))*(delta_t/60)*kwi(i)*exp(-kwo(i)*delta_t/60*(ttimes)) + est_fb(i)*imp) );
      
est_curves=est_curves((jitter+1):newRate:nTimes);   % downsampling step
     err = (tisscurve(1:nTimes,i) - est_curves);
      
     errsum = errsum + sum(err.*err);
     
     lambda=.1*errsum;
     lambda=0*errsum;   % changed to this Nov. 11 2005
     errsum = errsum + lambda*((t_delay(i)-estGlobalDelay)^2);



%clf 
%hold on
%plot(bldcurve,'r')
%plot(tmpbldcurve,'m')
%plot(est_curves(1:nTimes))
%%plot(tisscurve(:,i),'g')
%plot(tisscurve(:,i),'kx')
%drawnow
%t_delay

end

Fsum = errsum;
% disp('ERROR is '),Fsum
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


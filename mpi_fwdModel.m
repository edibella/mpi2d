%  mainly so can change how create fittedcurve  in only one place!!!

function [fittedcurve, tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, kwi, kwo, t_delay, est_fb);   % note all parameters are scalars!!!

nTimes=length(bldcurve);
tmp=0:nTimes-1;
ttimes = tmp';
imp=zeros(nTimes,1);
imp(1)=1.0;

% adding 7/2/04 to look at - even small say .05 fv cuts k1 in half..
%est_fb=0.07;
%t_delay=t_delay*0;

% t_delay just in the blood curve currently 4/04
% is t_delay always negative?  interp1 just doesn't give points < 0
    tmpbldcurve=bldcurve;
    tmpbldcurve(1-fix(t_delay):nTimes)=bldcurve(1:nTimes+fix(t_delay));
    tmpbldcurve(1:1-fix(t_delay))=0;
    t_delay=t_delay-fix(t_delay);
    tmpbldcurve = interp1(1:nTimes,tmpbldcurve,1+t_delay:nTimes+t_delay,'splines',0);
%    fittedcurve = conv(tmpbldcurve,((1-est_fb)*delta_t/60*kwi*exp(-kwo*delta_t/60*(ttimes)) + est_fb*imp) );
% changed 6/17/04 evrd to remove 1-est_fb: 
% and 10/17/05 to remove delta_t in front
% and change back 3/15/06!!!!
    fittedcurve = conv(tmpbldcurve,(delta_t/60*kwi*exp(-kwo*delta_t/60*(ttimes)) + est_fb*imp) );
% changed 9/13/05 to consider est_fb just models spillover. So no delay in it!
%    fittedcurve = conv(tmpbldcurve,(delta_t/60*kwi*exp(-kwo*delta_t/60*(ttimes)) ) );
%    fittedcurve=fittedcurve(1:nTimes) + est_fb*bldcurve;
    fittedcurve=fittedcurve(1:nTimes);


return;

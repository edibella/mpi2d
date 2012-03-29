function F=myfun_T1_subset_123_eq(x,xdata)
% This function is called by e.g. calT1_subset_P100609_eq.m
% and does not include the effects of T2*

N=24; % this must be correctly specified for the number of rays acquired in the sequence
TR=2.7;
TD=12;


numFrames=length(xdata(:,1));
total_err=0;

for ii=1:numFrames

   

  F(ii,:)=x(numFrames+2)*(((1-exp(-TR/x(ii)))/(1-cos(x(numFrames+1))*exp(-TR/x(ii))))+ ...
            (((cos(x(numFrames+1))*exp(-TR/x(ii))).^(xdata(ii,:)-1))*(1-(cos(x(numFrames+1))* ...
            exp(-TR/x(ii))).^N)/(N*(1-cos(x(numFrames+1))*exp(-TR/x(ii)))))* ...
            ((1-exp(-TD/x(ii)))-(1-exp(-TR/x(ii)))/(1-cos(x(numFrames+1))*exp(-TR/x(ii)))));

        
end


end



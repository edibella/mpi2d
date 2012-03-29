function [Gdconc Params]=calc_T1_subset_123_72ray(threeCurves,slice,precontrastT1)
 %%% ie. calc_T1_subset_123_eq('P100609',51,3,1950), or last entry can be left blank for series in which there was no contrast already "on board".
  
   
uload = threeCurves;

nonlcon=[];

s1(:,:)=uload(1:end,:); % for seperate frame
[frame dy]=size(s1(:,1));

N=24;   % rays per subset
TR=2.7;
TE=1.48;
TD=12;
%flipAngle=14;
   
    xdata1 = [1 25 49];
    xdata=repmat(xdata1,length(s1(:,1)),1);
    ydata = s1;
    
    %%%%% To separate the signal curves into 2 levels: low SI and high SI.
    %%%%% First solve fA and M0 for low SI curves to solve for T1(t), then use fixed fA and
    %%%%% bounds of M0 to estimate M0(t) and T1(t) for high SI regions.  In
    %%%%% the latter case, M0(t) varies likely due to T2* effects.
%     thresh=1.1; % choosing a reasonable threshold that is 10% higher than the tail of the highest subset s1 curve
%     m=(find(ydata(:,2)<thresh*mean(ydata(end-5:end,2))));
%     n=(find(ydata(:,2)>=thresh*mean(ydata(end-5:end,2))));
%     ydata1=ydata(m,:); % Subset s1 curves below the threshold
%     xdata1=xdata(m,:);
%     ydata2=ydata(n,:); % Subset s1 curves above the threshold
%     xdata2=xdata(n,:);
%     
    %s1=ydata1;
    
    xdata1=xdata;
    ydata1=ydata;
    
    
    
    [maxVal maxPos]=max(ydata1(:,1));
    [frame1 dy1]=size(ydata1(:,1));


    guess1=8.5*pi/180; % Initial Guess for the flip angle
    guess2=max(abs(s1(:,2))); % Initial Guess for M0
    x01 = [1500];  % Initial Guess for T1
    x0 = repmat(x01,length(ydata1(:,1)),1); % full vector of T1 guesses
    x0(length(ydata1(:,1))+1,1)=guess1; % concatenate the flip angle guess on the end
    x0(length(ydata1(:,1))+2,1)=guess2; % concatenate the M0 guess on the end
    
    LB(1:length(ydata1(:,1)),1)=0; % bounds for T1
    UB(1:length(ydata1(:,1)),1)=3000;
    LB(length(ydata1(:,1))+1,1)=5*pi/180; % bounds for flip angle
    UB(length(ydata1(:,1))+1,1)=20*pi/180;
    LB(length(ydata1(:,1))+2,1)=0; % bounds for M0
    UB(length(ydata1(:,1))+2,1)=2000000;


     [x] = lsqcurvefit(@myfun_T1_subset_123_eq,x0,xdata1,ydata1,LB,UB); 
     T1=x(1:length(ydata1(:,1)));
     fA=x(frame1+1,1); % in radians (must muliply by 180/pi to get in degrees)
     M0=x(frame1+2,1);
     


%%% Plugging in the T1, M0, and flip angle values to compute the signal for each time frame    
for jj=1:frame1
cSI_1(jj,1)=M0*(((1-exp(-TR/T1(jj,1)))./(1-cos(fA)*exp(-TR/T1(jj,1))))+ ...
            (((cos(fA)*exp(-TR/T1(jj,1))).^(xdata(1,1)-1))*(1-(cos(fA)*exp(-TR/T1(jj,1))).^N)/(N*(1-cos(fA)*exp(-TR/T1(jj,1)))))* ...
            ((1-exp(-TD/T1(jj,1)))-(1-exp(-TR/T1(jj,1)))./(1-cos(fA)*exp(-TR/T1(jj,1)))));
cSI_1(jj,2)=M0*(((1-exp(-TR/T1(jj,1)))./(1-cos(fA)*exp(-TR/T1(jj,1))))+ ...
            (((cos(fA)*exp(-TR/T1(jj,1))).^(xdata(1,2)-1))*(1-(cos(fA)*exp(-TR/T1(jj,1))).^N)/(N*(1-cos(fA)*exp(-TR/T1(jj,1)))))* ...
            ((1-exp(-TD/T1(jj,1)))-(1-exp(-TR/T1(jj,1)))./(1-cos(fA)*exp(-TR/T1(jj,1)))));
cSI_1(jj,3)=M0*(((1-exp(-TR/T1(jj,1)))./(1-cos(fA)*exp(-TR/T1(jj,1))))+ ...
            (((cos(fA)*exp(-TR/T1(jj,1))).^(xdata(1,3)-1))*(1-(cos(fA)*exp(-TR/T1(jj,1))).^N)/(N*(1-cos(fA)*exp(-TR/T1(jj,1)))))* ...
            ((1-exp(-TD/T1(jj,1)))-(1-exp(-TR/T1(jj,1)))./(1-cos(fA)*exp(-TR/T1(jj,1)))));
end

err=sum(sum((ydata-cSI_1).^2)); % toatal Chi-square error from all 3 curves
Params.T1=T1;
Params.M0=M0;
Params.fA=fA*180.0/pi;
Params.ChiSq=err;

max(T1);min(T1);

value(:,1)=[1:length(ydata1(:,1))]';
value(:,2)=T1;
R=0.0055; %L/mmol-ms 4 L/mmol-s
mean(T1(1:7));
if exist ('precontrastT1')
    Params.PrecontrastT1=precontrastT1;
    Gdconc(:,1)=((1./T1(:))-(1/precontrastT1))*(1/R);
else
   Gdconc(:,1)=((1./T1(:))-(1/mean(T1(1:7))))*(1/R);
   EstPrecontrast_T10=mean(T1(1:7))
   Params.EstPrecontrast=EstPrecontrast_T10;
end
value(:,3)=Gdconc;

figure; 
% subplot(1,3,1),plot(T1),axis square
% ylabel('T1 (ms)'),xlabel('frames')

%subplot(1,3,2)
hold on,plot(ydata1,'*'),plot(cSI_1),axis square
legend('S1(meas)','S2(meas)','S3(meas)','S1(est)','S2(est)','S3(est)')
ylabel('AIF subset curves w/ fits'),xlabel('frames')
%keyboard
% subplot(1,3,3), plot(1:frame,Gdconc,'ro-'),axis square
% ylabel('Gd (mmol/L)'); xlabel('frames');
% legend ('Gdconc') %,'Gdconc*10')

return
% AIF estimaion functions
% 8/31/10 - JUF

% input arguments 
%   data - [x,y,t] matrix
%   contours - 1X3 cell array with 
%       contours{1}=epi
%       contours{2}=endo
%       contours{3}=blood
%   time - 1XN array with timestamps 0:length(acquisition)

% output arguments
%   estAIF - 1XN array containing estimated AIF
%       other standard AMM output arguments may be added as needed

function estAif = tryMay2011CAMM(data,contours,time)

global satThresh

% set up CAMM parameters 
nCurves = 10;
nIter = 50;
useMeasAIF=1;
baselineScans=8;

settings.showPool = 0;
settings.showRunningAIF = 1;
settings.wordy = 0;
settings.plotAif = 1;
settings.plotErr = 1;
settings.findBestNumClus = 0;
settings.manSatLevel = 0;

% collect data

% convert signal to \Delta_SI
chis = zeros(size(data),'double');
S0 = squeeze(mean(single(data(:,:,1:baselineScans)),3));
for in=1:size(data,3)
    S = single(squeeze(data(:,:,in)));
    
    chis(:,:,in) = (S-S0);
end;

% collect tissue curves
sz=size(chis);
coords.endo=contours{2};
coords.epi=contours{1};
myoMask=getMask(chis,coords,1);

ctPool=zeros(sum(flatten(myoMask)),sz(3));
n=1;
for xi=1:sz(1)
    for yi=1:sz(2)
        if(myoMask(xi,yi))
            ctPool(n,:)=squeeze(chis(xi,yi,:));
            n=n+1;
        end
    end
end

keyboard

% collect AIF
x1poly=contours{3}(:,1); y1poly=contours{3}(:,2);
imagesc(squeeze(chis(:,:,fix(0.3*sz(3)))))
colormap gray
aifMask=roipoly(squeeze(chis(:,:,fix(0.3*sz(3)))),x1poly,y1poly);
h1=line(x1poly,y1poly);
set(h1,'Linewidth',2.5)
set(h1,'Color',[1 1 0]);

aifDat=zeros(sum(flatten(aifMask)),sz(3));
ii=1;
for xi=1:sz(1)
    for yi=1:sz(2)
        if(aifMask(xi,yi))
            aifDat(ii,:)=squeeze(chis(xi,yi,:))';
            ii=ii+1;
        end
    end
end

aif=median(aifDat);

% correct and load time
if(time(1)>0)
    time=time-time(1);
end
if(time(end)>5)
    time=time./60;
end
t=time(1:sz(3));

if (settings.showPool)
    if (sz(1) >= 100)
        N=100;
    else
        N=sz(1);
    end
    times=find(t > 0.05 & t < 5);
    figure(1)
    clf;
    hold on
    for ii=1:N
        plot(t(times),squeeze(ctPool(ii,times))+(ii-1)*0.3)
    end
    hold off
    pause(0.1)
end

if (settings.wordy)
    disp('Clustering ...')
end

if sz(1) == nCurves
    cts = ctPool;
elseif (settings.findBestNumClus)
    [~, cts]=calcBestNcurves(ctPool,nCurves,t,tInj);
else
    [~, cts]=kclusterSim(ctPool,nCurves);
end

% order of magnitude correction for scale
ctmax=max(flatten(cts));
if ctmax < 10
    ctscale=50/ctmax;
else
    ctscale=1;
end

cts=cts.*ctscale;
aif=aif.*ctscale;

tin=find(aif>10,1,'first');
tInj=t(tin-1);

% saturation selection
if (settings.manSatLevel)
    figure(1)
    plot(t,aif)
    satThresh=input('Saturation Level: ');
    mlim=max(aif);
else
    [mlim, mxind]=max(aif);
    [~, minind]=min(aif(mxind:mxind+30));
    minind=minind+mxind;
    [satThresh, satind]=max(aif(minind:end));
    satThresh=satThresh*1.4; %1.15;   % was 1.45 in 9/14 version, think Jacob changed 9/15/10
    satind=satind+minind;
    if (satind >= length(aif))
        disp('Saturation Level set to tail')
    end
end
    
meas_aif = aif;  % EVRD 9/20/10, so can plot out measured AIF for comparison
aif(aif >= satThresh)=satThresh;
cts=[aif;cts];
%

if (settings.showPool)
    figure(2)
    clf;
    plot(cts')
    pause(0.1)
end

% set initial guess
lb=[0 1 0 0 0 0 0 0 0 0 mlim];

options=optimset('TolFun',1e-8,'TolX',1e-8,'MaxIter',1e5,'Display','off');

popAifPar = [tInj, 4.051, 0.02766, 75.3613, 0.06213, 0.1331, 33.727, 0.5641,249.43,0.3546, 556.44];

if(size(t,1)>1)
    t=t';
end

if (useMeasAIF)
    measuredCppar = lsqcurvefit('satBlood',popAifPar,t,aif,lb,1000*ones(1,11),options);
    faifp = measuredCppar;
    avg=mean(aif(end-3:end));
else
    faifp = popAifPar;
    temp_aif=cpt(faifp,t);
    avg=mean(temp_aif(end-3:end));
end;

faifp(11)=500;     % was 550;  think Jacob changed 9/15/10
%faifp(11)=200;     % was 550;  think Jacob changed 9/15/10
faifp(11)=400

if (settings.wordy)
    disp('Estimating AIF ...')
end

% interpolation
% 10/4/10 - use median delta_t for interpolation                   
delta_t=t(2:end)-t(1:end-1);
temporalRes=median(delta_t);
tSamp=0:temporalRes:t(end);
ctsnew=zeros(size(cts,1),length(tSamp));
 for ii=1:size(cts,1)    
  ctsnew(ii,:)=interp1(t,cts(ii,:),tSamp,'cubic','extrap');  
 end

cts=ctsnew;     % end interpolation        
 
 

% estimation
[aifp,~,err] = SBCalternatingBlindEstimation(tSamp,cts,faifp,avg,nIter,settings);

disp('Fionished blind ')

if (settings.plotAif)
    figure(3); clf
    plot(t,aif,'r-',...
        t,cpt(aifp,t),'b.-',...
        t,(cpt(aifp,t)-aif),'go')
    hold on
    plot(t,meas_aif,'mx')
    legend('AIF measured, given', 'AIF estimated','Difference','satBldMeas')
end
if (settings.plotErr)
    figure(43)
    plot(err)
end
figure

estAif=cpt(aifp,t);
estAif=estAif./ctscale;

save tmpCAMMresults
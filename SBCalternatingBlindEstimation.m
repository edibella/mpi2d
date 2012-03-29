function [bpBest,PKpar,chis,runningAIFpar,runningPKpar] = SBCalternatingBlindEstimation(tmeas,Ctmeas,Cppar,avg,nIter,settings)

global PKpar satThresh

PKpar = [];

% set up lsqcurvefit params
lb=zeros(1,11);
lb(2)=1;
ub=5000*ones(1,7);
options=optimset('TolFun',1e-4,'TolX',1e-5,'MaxIter',1e4,'Display','off');

% initialize AIF
Cpguess = cpt(Cppar,tmeas);
 
Ncurves = size(Ctmeas,1);

% initialize kinetic parameters
PKpar(1,:)=[0 0 0 0];
for i=2:Ncurves
%     PKpar(i,:) = muraseFit(tmeas,Cpguess,Ctmeas(i,:),2);
%     PKpar(i,:) = muraseFit(tmeas,Cpguess,Ctmeas(i,:));
    PKpar(i,:) = lsqcurvefit('sbc_ct',[.8 3.7 .1 .05],[Cpguess tmeas],Ctmeas(i,:),[0 0 0 0],[4 10 1 2],options);
end;

chiCh=1;
chiSqold=1;

bpBest=Cppar;
errBest=1e6;
iterBest=1;
exectime = 0;

runningAIFpar=zeros(nIter,length(Cppar));
tic;
 
for i=1:nIter
    elapsedtimeOld=toc;
            
    if (settings.showRunningAIF)
        figure(2)
        clf;
        hold on;
        plot(tmeas,Ctmeas(1,:),'b.',tmeas,cpt(Cppar,tmeas),'r',tmeas,satBlood(Cppar,tmeas),'k');
        for ii=2:Ncurves
            plot(tmeas,Ctmeas(ii,:),'b.',tmeas,cpt(Cppar,tmeas),'r',tmeas,sbc_ct(PKpar(ii,:),[Cpguess,tmeas]),'k');
        end;
        hold off;
        title(['Iteration : ' num2str(i) ]);
        pause(.01);
%         pause
    end
    
    Cppar = lsqcurvefit(@SBCmodelCtPkFixed,Cppar,tmeas,Ctmeas,lb,ub,options);
    runningAIFpar(i,:)=Cppar;

    Cpguess = cpt(Cppar,tmeas);
    
    myavg=mean(Cpguess(end-3:end));
    scale=avg/myavg;
    Cpguess=Cpguess*scale;

    if (settings.showRunningAIF)
        figure(2)
        clf;
        hold on;
        plot(tmeas,Ctmeas(1,:),'b.',tmeas,Cpguess,'r',tmeas,satBlood(Cppar,tmeas),'k');
        for ii=2:Ncurves
            plot(tmeas,Ctmeas(ii,:),'b.',tmeas,Cpguess,'r',tmeas,sbc_ct(PKpar(ii,:),[Cpguess,tmeas]),'k');
        end;
        hold off;
        title(['Iteration : ' num2str(i) ]);
        pause(.01);
%         pause
    end

 
    for ii=2:Ncurves
%         PKpar(i,:) = muraseFit(tmeas,Cpguess,Ctmeas(i,:),2);
%         PKpar(i,:) = muraseFit(tmeas,Cpguess,Ctmeas(i,:));
        PKpar(ii,:) = lsqcurvefit('sbc_ct',PKpar(ii,:),[Cpguess tmeas],Ctmeas(ii,:),[0 0 0 0],[5 25 1 2],options);
    end;
    
    runningPKpar(i,:,:)=PKpar;
    
    for ii=2:Ncurves
        CtEst(ii,:)=sbc_ct(PKpar(ii,:),[Cpguess,tmeas]);
        chi(ii)=sum((CtEst(ii,:)-Ctmeas(ii,:)).^2);
    end
        
    chiSq=sum(chi);
    
    chiCh=abs(chiSqold-chiSq);
    
    chiSqold=chiSq;
    chis(i)=chiSq;
    
    % update 'best' values
    if chis(i) <= errBest
        errBest=chis(i);
        iterBest=i;
        bpBest=Cppar;
    end
    
    % reset if off course
    if chis(i) >= 2* errBest
        Cppar = bpBest+(0.1*rand(size(bpBest))-0.05).*bpBest;
        PKpar = squeeze(runningPKpar(iterBest,:,:));
        disp(['Resetting from iteration: ',num2str(i),' to iteration: ',num2str(iterBest)]);
    end
%           
%     if chiCh <= 0.001*chis(i) &&  chiSqold <= errBest*1.001 && i >= floor(nIter/4)
%         break
%     end
% EVRD 12/20/11, see if helps to not stop so early
    
    if chiCh <= 0.0001*chis(i) &&  chiSqold <= errBest*1.001 && i >= floor(nIter/4)
        break
    end
    
    elapsedtime = toc;
    
    exectime=elapsedtime-elapsedtimeOld;
    
    if(exectime >= 250)
        break
    end
    
    if (settings.wordy)
        disp(['Iteration: ',num2str(i),'; Execution time: ', num2str(exectime)]);
    end
        
end;

return;

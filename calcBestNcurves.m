function [indx cts]=calcBestNcurves(ctPool,nCurves,t,tInj)

aif = cpt([tInj, 4.051, 0.02766, 75.3613, 0.06213, 0.1331, 33.727, 0.5641,249.43,0.3546, 556.44],t);

li=1;
N=4:4:nCurves*4;
sd=zeros(N,1);
s=sd;
d=sd;
ncur=sd;
for ii=6:4:nCurves*4
    
    snrDiv = zeros(10,1);
    div = zeros(10,1);
    snr = zeros(10,1);

    disp(['      Testing ',num2str(ii),' clusters...']);
    for zi=1:10
        
        [~, tmpcts] = kclusterSim(ctPool,ii);
        
        tcMax=zeros(1,ii);
        noise=zeros(1,ii);
        
        tmp2=zeros(ii,4);
        nC = size(tmpcts,1);
        for ki=1:nC
            tmp2(ki,:) = muraseFit(t,aif,tmpcts(ki,:));
            tmp2 = max(tmp2,0);
            measKps(ki,:) = tmp2(ki,1:3) ./ [1.2 1 .2];
        end
        if(max(measKps(:,3))==0 | max(measKps(:,2))==0)
            A=0;
        else
            [~, A]=convhull(measKps(:,2),measKps(:,3));
        end
%         div(zi)=3*A;
        div(zi)=1.421*A;
        % plot(NmeasKps(K,2),NmeasKps(K,3),'r-',NmeasKps(:,2),NmeasKps(:,3),'b+')
        for ji=1:nC
            tcMax(ji)=max(tmpcts(ji,:));
            noise(ji)=std(tmpcts(ji,1:5));
        end
        
        signal=nanquantile(tcMax',.75);
        
        snr(zi)=.0608*(signal/mean(noise));
%         snr(zi)=.03*(signal/mean(noise));
        
        snrDiv(zi)=snr(zi)+div(zi);
        
    end
    ncur(li)=ii;
    s(li)=mean(snr);
    d(li)=mean(div);
    sd(li)=mean(snrDiv);
    li=li+1;
    
end

bestNc = ncur(sd==max(sd));
if length(bestNc > 1)
    bestNc=bestNc(end);
end

disp(['     Optimal number of clusters: ',num2str(bestNc)]);
plot(sd,'b')
hold on
plot(s,'g')
plot(d,'r')
hold off
pause(1)
[indx cts]=kclusterSim(ctPool,bestNc);


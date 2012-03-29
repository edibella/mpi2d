%%%%%
%
% Cluster
%
%%%%%

% input
% ctvals - Tissue Curve Matrix 
% numclu - number of clusters

% output
% indx   - Cluster index
% cts2   - Averaged tissue curves from each cluster

%function [indx cts2 cts1]=kclusterSim(ctvals,numclu)   

function [indx cts2]=kclusterSim(ctvals,numclu)   

[nCurves lent]=size(ctvals);

%[indx c2]=k_means(ctvals,numclu);
[indx c2]=kmeans(ctvals,numclu,'emptyaction','singleton');

inx=[];
for ii=1:numclu
	if (max(indx==ii)==1)
		inx=[inx,ii];
	end
end

cts2=zeros(length(inx),lent);


for ii=1:length(inx)
	ji=inx(ii);
    if sum(indx==ji) > 1
        cts2(ji,:)=mean(ctvals(indx==ji,:));
    else
        cts2(ji,:)=ctvals(indx==ji,:);
    end
% 	figure(1)
% 	clf;
% 	plot(transpose(kpMat(indx==ji,:)))
		
	figure(2)
	clf;
	plot(transpose(ctvals(indx==ji,:)))
	hold on
	plot(cts2(ji,:),'+')
% 	hold off
% 	title(num2str(ji))
% 	pause(0.5)
% 	
% 	noise(ji)=std(cts2(ji,1:15));
% 	avgsnr=mean(std(ctvals(indx==ji,1:15)));
% 	
% 	fprintf('Average Data Noise: %5.4f, Clustered Noise %5.4f, Improvement %5.4f \n',avgsnr,noise(ji),avgsnr-noise(ji));
% 	fprintf('Ncurves: %g, Factor %2.1f \n',sum(indx==ji),avgsnr/noise(ji));
end

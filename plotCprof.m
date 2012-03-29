function plotCprof(vector, studyNum, vectorDescriptor, studyDescriptor, lineType, colorVector)

%clf
set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]); hold on
xlabel(strcat(vectorDescriptor,' Region Number'));
ylabel('Flow value')

for ii=1:size(vector,2)
      plot(vector(:,ii), strcat(lineType,colorVector(ii)),'linewidth',2);
end

[numAzimuthalRegions numSlices]=size(vector);

avgkwi=mean(mean(vector));
% maybe better, compute avg of highest half:
for islice=0:numSlices-1
  for jj=1:numAzimuthalRegions
     aa(jj+islice*numAzimuthalRegions)=vector(jj,islice+1);
  end
end
  bb=sort(aa);
  avgOfVector=mean(bb(round(0.5*length(bb)):length(bb)));
%  plot(avgOfVector*ones(1,numAzimuthalRegions),':','linewidth',1);
  Ylim=get(gca,'YLim');
  Ylim(1)=0.0;
  Xlim=get(gca,'Xlim')
%  axis([Xlim Ylim]);

%legend('Slice 1 - most apical','Slice 2', 'Slice 3 - most basal', 'avg of upper half')
if size(vector,2)==3
    legend('Slice 1 - most apical','Slice 2', 'Slice 3 - most basal')
elseif size(vector,2)==4
    legend('Slice 1 - most apical','Slice 2', 'Slice 3', 'Slice 4 - most basal')
else  %if size(vector,2)==5
    legend('Slice 1 - most apical','Slice 2', 'Slice 3', 'Slice 4', 'Slice 5 - most basal')
end

disp('Note these slice numbers do not correspond to acquisition!!!')
title(strcat(studyDescriptor,' Series ',int2str(studyNum)))

drawnow


return;

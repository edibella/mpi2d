function showcurves(curves, yaxislabel)

[nRegs, nTimes] = size(curves);
bldcurve=curves(1,:);
tisscurves=curves(2:nRegs,:);
tisscurve=tisscurves';
%nTimes=length(bldcurve);

clf
hold on

set(gca,'FontSize',16)
set(gcf,'Color',[1 1 1])   % to get rid of grey border
%subplot(2,1,1)
hold on
set(gca,'FontSize',16)
%plot(1:nTimes,bldcurve,'-r','Linewidth',2)


if nRegs==17   % assume this is endo and epi for now 
   plot(1:nTimes,tisscurve(:,1:floor(nRegs/2)),'Linewidth',2)
%axis([0 30 -0.15 0.85])
   Ylim=get(gca,'YLim');
   Xlim=get(gca,'Xlim')
   axis([Xlim Ylim]);
%yt = get(ax,'YTick');  % Y-tick locations
%xc = get(ax,'XColor'); % Color of X-axis
%set(gca,'XTickMode','manual', ...  % Fix the tick limits and
%        'XLimMode','manual', ...   % locations
set(gca,'XTickLabel',[' ']);
%legend('Blood ROI', 'Septal','Anterior','Lateral','Inferior')
legend('Reg 1','Reg 2','Reg 3','Reg 4','5','6','7','8')
title('Endocardial curves')
figure; set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
   plot(1:nTimes,tisscurve(:,floor(nRegs/2)+1:nRegs-1),'Linewidth',2)
   legend('Reg 1','Reg 2','Reg 3','Reg 4','5','6','7','8')
   title('Epicardial curves')
   axis([Xlim Ylim]);
figure; set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
end   %if nRegs=17


%subplot(2,1,2)
clf
hold on
set(gca,'FontSize',16)
%plot(1:nTimes,bldcurve,'-xr','Linewidth',2)
%plot(1:nTimes,tisscurve(:,floor(nRegs/2)+1:nRegs-1),'x','Linewidth',0.5)
%plot(1:nTimes,tisscurve(:,floor(nRegs/2)+1:nRegs-1),'Linewidth',2)
plot(1:nTimes,tisscurve,'Linewidth',2)
plot(1:nTimes,tisscurve,'o','Linewidth',1)
%axis([0 30 -0.15 0.85])
ylabel(yaxislabel)
xlabel('Frame number')

%legend('Blood', 'Reg 5','Reg 6','Reg 7','Reg 8')
%legend('Blood', '1','2','3','4','5','6','7','8')

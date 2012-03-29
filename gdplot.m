function scaleupGd=gdplot(Gdcon,deltaSIcurves)


if length(deltaSIcurves(1,:))>1
    deltaSI=deltaSIcurves(1,:);
else
    deltaSI=deltaSIcurves;
end

[dt]=length(Gdcon);
[dt2]=length(deltaSI);

if dt> dt2
    timeframe=dt2;
else
    timeframe=dt;
end

scaleb=timeframe;
%===============
%the point to scale up to deltaSIcurve
%scaleb=32;
%===============
Gd=Gdcon-mean(Gdcon(1:5));
scale2=mean(deltaSI(scaleb-5:scaleb-2))/mean(Gd(scaleb-5:scaleb-2));
scaleupGd=Gd*scale2;
figure(21),hold on, plot(1:timeframe,deltaSI(1:timeframe),'b-*')
plot(1:timeframe,Gd*scale2,'g-o')

%%==========================

legend ('deltaSI','Gdcon')

xlabel('Image Frames'),ylabel('deltaSI')
%title('P052709')


end % function
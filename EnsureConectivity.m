function connectedMyo = EnsureConectivity(myo,Center)
    [sx sy] = size(myo);
    %now make contours that extracts out the myocardium
    [r,c] = find(myo>0);
    [theta,radius] = cart2pol(r-Center(1),c- Center(2));
    if(isempty(theta))
        connectedMyo = [];
        return;
    end
    theta = mod(theta,2*pi);
    thetavals = 0:.04:(2*pi);
    SampledRs = [];
    thetas = [];
    dt = 0.1;
    for thetai=1:length(thetavals)
        thetaval = thetavals(thetai);
        meanR = mean(radius((theta< thetaval+dt).*(theta> thetaval-dt) > 0));
        if(~(isempty(meanR) || isnan(meanR)))
            SampledRs = [SampledRs meanR];
            thetas = [thetas thetaval];
        end
    end
    SampledRs = horzcat(SampledRs,SampledRs,SampledRs);
    thetas = [thetas (thetas+2*pi) (thetas + 4*pi)];
    thetavals = 0:.01:(2*pi);
    desiredThetas = thetavals;
    desiredThetas = [desiredThetas (desiredThetas+2*pi) (desiredThetas + 4*pi)];
    interpolatedRs = interp1(thetas,SampledRs,desiredThetas,'pchip');
    interpolatedRs = interpolatedRs((length(thetavals)+1):(2*length(thetavals)));
    connectedMyo = myo;
    for thetai = 1:length(thetavals)
        thetaval = thetavals(thetai);
        meanR = interpolatedRs(thetai);
        x = round(meanR*cos(thetaval) + Center(1));
        y = round(meanR*sin(thetaval) + Center(2));
        y = max(min(y,sy-1),1);
        x = max(min(x,sx-1),1);
        connectedMyo((x):(x+1),(y):(y+1)) = 1;
        %figure(1),clf, imagesc(connectedMyo), hold on, plot(y,x,'go'),pause(.05);
    end
    
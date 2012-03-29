function [newMask_out] = maintainMinumumConnectivity(newMask)
newMask_out = newMask;
[r,c] = find(newMask>0);
middle = mean([r c]);
[mytheta, myradius] = cart2pol(r-middle(1), c-middle(2));
mytheta = [mytheta;(mytheta+2*pi)];
myradius = [myradius;myradius];
fixedSomething = 0;
dtheta = .05;
runningRadius = 0;
for theta = (-pi):dtheta:(4*pi)
    lowerTheta = theta-dtheta;
    upperTheta = theta+dtheta;
    indicies = (mytheta>= lowerTheta) .*(mytheta<= upperTheta);
    if(sum(indicies) == 0)
        fixedSomething = 1;
        for theta2 = lowerTheta:upperTheta
            mytheta = [mytheta; theta2];
            myradius = [myradius; runningRadius];
        end
    else
        runningRadius = mean(myradius(indicies>0));
        mythinradius(indicies>0) = runningRadius;
    end
end

if(fixedSomething)
    [x,y] = pol2cart(mytheta,myradius);
    x = round(x + middle(1));
    y = round(y + middle(2));
    minimumConnectivity_out = zeros(size(minimumConnectivity));
    newMask_out = newMask;
    for i=1:length(x)
        if(newMask_out(x(i),y(i)) ~= 1)
            newMask_out(x(i),y(i)) = 1;
        end
    end
end
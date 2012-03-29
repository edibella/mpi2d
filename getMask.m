function myMask=getMask(data,coords,st_ang)

[lenx leny lent]=size(data);

%% Draw end and epi
% figure

tmpn=coords.endo;
tmpp=coords.epi;
imagesc(squeeze(data(:,:,fix(0.3*lent))))
colormap gray

x1poly=tmpn(:,1); y1poly=tmpn(:,2);
bw1=roipoly(squeeze(data(:,:,fix(0.3*lent))),x1poly,y1poly);
h1=line(x1poly,y1poly);
set(h1,'Color',[0 1 0]); 
set(h1,'Linewidth',2.5)
%set(h1,'Color',[1 1 0]); set(h1,'Linestyle','x'); set (h1,'Linestyle','-'); 
set(h1,'Color',[1 1 0]);
[I, J]=find(bw1);
center1=[mean(x1poly) mean(y1poly)];


x2poly=tmpp(:,1); y2poly=tmpp(:,2);
bw2=roipoly(squeeze(data(:,:,fix(0.3*lent))),x2poly,y2poly);
h2=line(x2poly,y2poly);
set(h2,'Linewidth',2.5)
center2=[mean(x2poly) mean(y2poly)];

start_angle=st_ang;
while start_angle >360
    start_angle=start_angle-360;
end
    

%% Draw regions... create region wide-mask

center=(center1+center2)/2;
xCenter = center(1);
yCenter = center(2);

myo = double(bw2)-double(bw1);   % check all are 0 or 1 ... (2nd mask encompasses first!)
[Y, X] = find(myo);
nX     = length(X);

% numAzimuthalRegions = input('Number of regions: ');
numAzimuthalRegions = 1;

dAng  = 360 / numAzimuthalRegions;
mask = zeros(lenx, leny, numAzimuthalRegions);
Angles =  atan2(-(Y - yCenter), X - xCenter) * 180/pi + 180;
for i = 1 : numAzimuthalRegions
    clear RoiX RoiY;
    cn = 0;
    Ang0 = (i - 1) * dAng + start_angle;
    Ang1 = i * dAng       + start_angle;
    for j = 1 : nX
        if( ((Angles(j) > Ang0) && (Angles(j) <= Ang1)) || ((Angles(j) + 360 > Ang0) && (Angles(j) + 360 <= Ang1)))
            cn = cn + 1;
            RoiX(cn) = X(j);
            RoiY(cn) = Y(j);
        end
    end

    for n = 1 : cn
        mask(RoiY(n), RoiX(n), i) = 1;

    end

    %   plot(RoiX, RoiY, '.', 'Color', [0.5, 0.7, 0.9]);
    rr = lenx/2;
    xx = [xCenter,    -rr * cos(Ang1*pi/180) + xCenter] ;
    yy = [yCenter,     rr * sin(Ang1*pi/180) + yCenter] ;

    % find intersection with contour:


    tmpp=-atan2(y2poly-yCenter,x2poly-xCenter)*180/pi+180;

    [closestAngle, indexIntersect]=sort(abs(tmpp-Ang1));
    radDistToEpi=sqrt((x2poly(indexIntersect(1))-xCenter)*(x2poly(indexIntersect(1))-xCenter)+(y2poly(indexIntersect(1))-yCenter)*(y2poly(indexIntersect(1))-yCenter));

    xx=[xCenter, xCenter-1.4*(radDistToEpi*cos(Ang1*pi/180)) ] ;
    yy = [yCenter, yCenter+1.4*(radDistToEpi*sin(Ang1*pi/180))];
    hh=line(xx,yy);

    %        find xx in myo
    set(hh,'Color',[0 0 1])
    set(hh,'Linewidth',3)

    xx=[xCenter-1.4*(radDistToEpi*cos(Ang1*pi/180+pi/numAzimuthalRegions)) ] ;
    yy = [yCenter+1.4*(radDistToEpi*sin(Ang1*pi/180+pi/numAzimuthalRegions))];
    if i==numAzimuthalRegions
        labelString=sprintf('%d',1);
    else
        labelString=sprintf('%d',i+1);
    end
    ht=text(xx,yy,labelString);
    set(ht,'FontSize',18,'Color',[1 1 1])

    t = mask(:,:,i);
    N = sum(sum(t));
    pause(2)
    
end
myMask(:,:)=mask;


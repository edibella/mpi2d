function newContour = morph(oldcontour, radius)

%find limits of contour given
maxx = max(oldcontour(1,:));
maxy = max(oldcontour(2,:));
space = zeros(maxx+2*radius, maxy + 2*radius);

%find a mask enclosing the specified contour
[BW x y] = roipoly(space,oldcontour(1,:),oldcontour(2,:));

%perform a morph based on it
MorphType = {'dilate','erode'};
MorphType = MorphType{(radius>0)+1};
BW = bwmorph(v.epi.mask,MorphType,abs(radius));

%find a contour based on the new contour
[r c] = find(v.endo.mask);
[minr minri] = min(r);
v.endo.boundary = bwtraceboundary(BW,[r(minri) c(minri)],'N');
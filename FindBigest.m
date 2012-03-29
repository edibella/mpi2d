function finalMask = FindBigest(segmented)
parts = max(segmented(:));
se = strel('diamond',1);
newmaskValue = parts;
receptacle = zeros(size(segmented));
sizes = [];
for i=1:parts
    mymask = (segmented==i);
    while(sum(mymask(:)))
        newmaskValue = newmaskValue + 1;
        [r, c] = find(mymask);
        BW = imfill(mymask,[r(1) c(1)]);
        sizes(newmaskValue) = sum(BW(:));
        receptacle = receptacle + BW*newmaskValue;
        mymask = mymask - BW;
    end
end
finalMask = zeros(size(mymask));
for i=1:4
    [bestValue imax] = max(sizes);
    sizes(imax) = 0;
    finalMask = finalMask + i*(receptacle==imax);
end
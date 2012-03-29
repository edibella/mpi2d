function smoothedImg = smoothimg(paddedmovie, paddedoffset, index, weightingmatrix)
durationofWeightingMatrix = size(weightingmatrix,3);
weightingoffset = round((durationofWeightingMatrix-1)/2);
smoothedImg = sum(paddedmovie(:,:,paddedoffset + index + (-weightingoffset:weightingoffset)).*weightingmatrix,3);
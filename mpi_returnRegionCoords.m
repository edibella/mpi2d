function [RoiXcoords, RoiYcoords]=mpi_regions(ireg, Xcoords, Ycoords, numAzimuthalRegions, start_angle, pixelSizeRatio, xCenter,yCenter)

   nX=length(Xcoords);
   clear RoiXcoords RoiYcoords;
   cn = 0;
   dAng  = 360 / numAzimuthalRegions;
   Ang0 = (ireg - 1) * dAng + start_angle;
   Ang1 = ireg * dAng       + start_angle;
   Angles =  atan2(-(Ycoords - yCenter), pixelSizeRatio*(Xcoords - xCenter)) * 180/pi + 180; 
   for j = 1 : nX
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 < Ang1)) | ((Angles(j)+720 > Ang0) & (Angles(j)+720 < Ang1)))
         cn = cn + 1;
         RoiXcoords(cn) = Xcoords(j);
         RoiYcoords(cn) = Ycoords(j);
      end
   end

return;

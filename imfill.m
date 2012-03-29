function BW = imfill(img, startLoc)
pointsToCheck = startLoc;
[sx sy] = size(img);
BW = zeros(size(img));
BW(startLoc(1), startLoc(2)) = 1;
while(length(pointsToCheck) > 0)
    point = pointsToCheck(1,:);
    pointsToCheck(1,:) = [];
    BW(point(1), point(2)) = 1;
    
    if((point(1) < sx) && (img(point(1)+1,point(2)) == 1) && (BW(point(1)+1,point(2)) ~= 1))
        BW(point(1)+1,point(2)) = 1;
        pointsToCheck = vertcat(pointsToCheck,[point(1)+1 point(2)]);
    end
    if((point(2) < sy) && (img(point(1),point(2)+1) == 1) && (BW(point(1),point(2)+1) ~= 1))
        BW(point(1),point(2)+1) = 1;
        pointsToCheck = vertcat(pointsToCheck,[point(1) point(2)+1]);
    end
    if((point(1) > 1) && (img(point(1)-1,point(2)) == 1) && (BW(point(1)-1,point(2)) ~= 1))
        BW(point(1)-1,point(2)) = 1;
        pointsToCheck = vertcat(pointsToCheck,[point(1)-1 point(2)]);
    end
    if((point(2) > 1) && (img(point(1),point(2)-1) == 1) && (BW(point(1),point(2)-1) ~= 1))
        BW(point(1),point(2)-1) = 1;
        pointsToCheck = vertcat(pointsToCheck,[point(1) point(2)-1]);
    end
    if((point(1) < sx) && (point(2) < sy) && (img(point(1)+1,point(2)+1) == 1) && (BW(point(1)+1,point(2)+1) ~= 1))
        BW(point(1)+1,point(2)+1) = 1;
        pointsToCheck = vertcat(pointsToCheck,[point(1)+1 point(2)+1]);
    end
    if((point(1) > 1) && (point(2) < sy) && (img(point(1)-1,point(2)+1) == 1) && (BW(point(1)-1,point(2)+1) ~= 1))
        BW(point(1)-1,point(2)+1) = 1;
        pointsToCheck = vertcat(pointsToCheck,[point(1)-1 point(2)+1]);
    end
    if((point(1) > 1) && (point(2) > 1) && (img(point(1)-1,point(2)-1) == 1) && (BW(point(1)-1,point(2)-1) ~= 1))
        BW(point(1)-1,point(2)-1) = 1;
        pointsToCheck = vertcat(pointsToCheck,[point(1)-1 point(2)-1]);
    end
    if((point(1) < sx) && (point(2) > 1) && (img(point(1)+1,point(2)-1) == 1) && (BW(point(1)+1,point(2)-1) ~= 1))
        BW(point(1)+1,point(2)-1) = 1;
        pointsToCheck = vertcat(pointsToCheck,[point(1)+1 point(2)-1]);
    end
end

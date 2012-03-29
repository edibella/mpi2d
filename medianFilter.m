function IM = medianFilter(img)
[sx sy] = size(img);
for i=2:(sx-1)
    for j=2:(sy-1)
        temp = img((i-1):(i+1),(j-1):(j+1));
        IM(i,j) = median(temp(:));
    end
end
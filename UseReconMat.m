cd('/v/raid1/bmatthew/MRIdata/Cardiac/Trio/P090309/Processing')
seriesNumbers = 22:24
for slice=1:length(seriesNumbers)
    ReconMat(51, seriesNumbers(slice),slice);
end
seriesNumbers = 27:29
for slice=1:length(seriesNumbers)
    ReconMat(55, seriesNumbers(slice),slice);
end

seriesNumbers = 31:32
for slice=1:length(seriesNumbers)
    ReconMat(57, seriesNumbers(slice),slice);
end

seriesNumbers = 64:65
for slice=1:length(seriesNumbers)
    ReconMat(76, seriesNumbers(slice),slice);
end

function Ktrans = getKtrans(f)
fid=fopen(f,'r');
if(fid < 0)
    disp(['Couldn''t open file: ' f]);
    Ktrans = -1;
    return
end
Ktrans = [];
tline = fgetl(fid);
while ischar(tline)
    if(~isempty(strfind(tline,'Ktrans  ')))
        A = sscanf(tline,'Ktrans  %d     %f');
        Ktrans(A(1)) = A(2);
    end
    if(~isempty(strfind(tline,'flow  ')))
        A = sscanf(tline,'flow  %d     %f');
        Ktrans(A(1)) = A(2);
    end
    tline = fgetl(fid);
end
fclose(fid);
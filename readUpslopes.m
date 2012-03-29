function [upslope]=readUpslopes(filename)

% last two may or may not be there...

fid=fopen(filename ,'r' );
a =fscanf(fid,'%s',3);
a =fscanf(fid,'%s',1);
delta_t=str2num(a(9:14));

indexUpslope=1;
while ~isempty(a)
  a =fscanf(fid,'%s',1);
  switch(a)
     case {'upslope'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
        aa=str2num(a);
        upslope(indexUpslope)=aa; indexUpslope=indexUpslope+1;
     otherwise
  end
end

return

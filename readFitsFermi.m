function [FLOW_0, Scalar_F,ve, est_fb, t_delay, Tzero, fval, delta_t, spillover]=readFitsFermi(filename)

% last two may or may not be there...

fid=fopen(filename ,'r' );
a =fscanf(fid,'%s',3);
a =fscanf(fid,'%s',1);
delta_t=str2num(a(9:14));

indexFLOW_0=1; indexScalar_F=1; indexve=1; indexfv=1; indexTzero=1; indexfval=1;
indext_delay=1; indexspillover=1; indexest_fb=1;
while ~isempty(a)
  a =fscanf(fid,'%s',1);
  switch(a)
     case {'FLOW_0'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
        aa=str2num(a);
        FLOW_0(indexFLOW_0)=aa; indexFLOW_0=indexFLOW_0+1;
     case {'Scalar_F'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
        aa=str2num(a);
        Scalar_F(indexScalar_F)=aa; indexScalar_F=indexScalar_F+1;
     case {'kwo'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        ve(indexve)=aa;  indexve=indexve+1;
     case {'fv'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        est_fb(indexfv)=aa;  indexfv=indexfv+1;
     case {'t0'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        t_delay(indext_delay)=aa;  indext_delay=indext_delay+1;
     case {'Tzero'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        Tzero(indexTzero)=aa;  indexTzero=indexTzero+1;
     case {'est_fb'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        est_fb(indexest_fb)=aa;  indexest_fb=indexest_fb+1;
      case {'spillover'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        spillover(indexspillover)=aa;  indexspillover=indexspillover+1;
     case {'fval'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        fval(indexfval)=aa;  indexfval=indexfval+1;
     otherwise
  end
end
fclose(fid);
% for global case
if ~exist('spillover')
  spillover=zeros(size(FLOW_0));
end
if ~exist('est_fb')
    est_fb=zeros(size(FLOW_0));
end


Hct=0.45;
%ve=kwo;  % not at all equivalent, but for convenience write out thsi way
%est_fb=Tzero;   % changed 4/06

%params=[flow,ve, est_fb, t_delay, fval, delta_t, spillover]
return

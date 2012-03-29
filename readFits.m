function [Ktrans,ve, est_fb, t_delay, fval, delta_t, spillover]=readFits(filename)

% last two may or may not be there...
fid=fopen(filename ,'r' );
if(fid < 0)
    x = 0;
end
a =fscanf(fid,'%s',3);
a =fscanf(fid,'%s',1);
delta_t=str2num(a(9:14));

indexflow=1; indexve=1; indexfv=1; indext0=1; indexfval=1;
indext_delayVp=1; indexspillover=1;
while ~isempty(a)
  a =fscanf(fid,'%s',1);
  switch(a)
     case {'Ktrans'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
        aa=str2num(a);
        Ktrans(indexflow)=aa; indexflow=indexflow+1;
     case {'ve'}
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
        t_delay(indext0)=aa;  indext0=indext0+1;
     case {'t_delayVp'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        t_delayVp(indext_delayVp)=aa;  indext_delayVp=indext_delayVp+1;
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
%%% Below was commented to call "Ktrans" from the output 'flowvalues...' file, instead of "flow"
% % % for global case
% % if ~exist('spillover')
% %   spillover=zeros(size(flow));
% % end
% % if ~exist('est_fb')
% %     est_fb=zeros(size(flow));
% % end
% % if ~exist('t_delay')
% %     t_delay=zeros(size(flow));
% % end

% for global case
if ~exist('spillover')
  spillover=zeros(size(Ktrans));
end
if ~exist('est_fb')
    est_fb=zeros(size(Ktrans));
end
if ~exist('t_delay')
    t_delay=zeros(size(Ktrans));
end
Hct=0.45;
%kwo=(1-Hct)*(flow*0.5)./ve;
kwo=(1-Hct)*(Ktrans)./ve;

%params=[flow ve est_fb t_delay fval];
return

function [flow,ve, est_fb, t_delay, fval, delta_t, spillover]=readFitsFile(filename,modelType)

if strcmp(modelType,'fermi')
   filename=strcat(filename,'.fermi.txt')
   [flow ve est_fb t_delay fval delta_t]=readFitsFermi(filename);
   disp('fermi......')
   return
elseif strcmp(modelType,'fermi_constrDelay')
   filename=strcat(filename,'.fermi.constrDelay.txt')
   [flow ve est_fb t_delay fval delta_t]=readFitsFermi(filename);
   return
elseif strcmp(modelType,'fermiFull')
   filename=strcat(filename,'.fermiFull.txt')
   [flow ve est_fb t_delay fval delta_t]=readFitsFermi(filename);
   return
elseif strcmp(modelType,'fullModel2')
   filename=strcat(filename,'.fullModel2.txt')
elseif strcmp(modelType,'globaldelay')
   filename=strcat(filename,'.globaldelay.txt')
elseif strcmp(modelType,'noBlood')
   filename=strcat(filename,'.noBlood.txt')
elseif strcmp(modelType,'noBlood1')
   filename=strcat(filename,'.noBlood1.txt')
elseif strcmp(modelType,'full')
   filename=strcat(filename,'.full.txt')
 elseif strcmp(modelType,'globaldelayAndVp')
    filename=strcat(filename,'.globaldelayAndVp.txt')
elseif strcmp(modelType,'funcfit')
   filename=strcat(filename,'.funcfit.txt')
elseif strcmp(modelType,'noBlood2')
   filename=strcat(filename,'.noBlood2.txt')
elseif strcmp(modelType,'fermi2')
   filename=strcat(filename,'.fermi2.txt')
   [flow ve est_fb t_delay fval delta_t]=readFitsFermi(filename)
   return;
elseif strcmp(modelType,'fermi_constrDelayBld')
   filename=strcat(filename,'.fermi.constrDelayBld.txt')
   [flow ve est_fb t_delay fval delta_t]=readFitsFermi(filename);
   return
elseif strcmp(modelType,'globaldelay2')
   filename=strcat(filename,'.globaldelay2.txt')
elseif strcmp(modelType,'globaldelayAndVp2')
   filename=strcat(filename,'.globaldelayAndVp2.txt')

elseif strcmp(modelType,'upslopes')
   filename=strcat(filename,'.upslopes.txt')
   [flow]=readUpslopes(filename);
   fval=zeros(1,nRegs);
   return
end
%else
%   [flow ve est_fb t_delay fval delta_t]=readFits(filename);
%end



% last two may or may not be there...
filename

fid=fopen(filename ,'r' );
a =fscanf(fid,'%s',3);
a =fscanf(fid,'%s',1);
delta_t=str2num(a(9:14));

indexflow=1; indexve=1; indexfv=1; indext0=1; indexfval=1;
indext_delayVp=1; indexspillover=1;
while ~isempty(a)
  a =fscanf(fid,'%s',1);
  switch(a)
     case {'flow'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
        aa=str2num(a);
        flow(indexflow)=aa; indexflow=indexflow+1;
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

% for global case
if ~exist('spillover')
  spillover=zeros(size(flow));
end
if ~exist('est_fb')
    est_fb=zeros(size(flow));
end
if ~exist('t_delay')
    t_delay=zeros(size(flow));
end


Hct=0.45;
kwo=(1-Hct)*(flow*0.5)./ve;


%params=[flow ve est_fb t_delay fval];
return

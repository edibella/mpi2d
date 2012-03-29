function printFitParams(xx,fval,sliceNum, studyNum, outpath, delta_t, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF,flagTimeStamps,modelType)

if nargin==0,
    error('arguments needed for printFitParams');
end
%if ~exist('numSkip'),   % if manually skipping an intital bump desired
%   numSkip=0;
%end

if(flagPixelwise);
    for ii=1:length(xx)
        kwi(ii) = xx(ii,1);
        kwo(ii)= xx(ii,2);
        %PartCoeff(ii)= xx(ii,2);% this was changed to to replace kwo above (NP 10/10/06)
        est_fb(ii)=fixedVp;
        if fixedDelay==99
            t_delay(ii) = xx(ii,3);
        else
            t_delay(ii)=fixedDelay;   % came from estGlobalDelay
            if fixedVp==99
                est_fb(ii)= xx(ii,3);
            end
        end
        if fixedVp==99 & fixedDelay==99
            est_fb(ii)= xx(ii,4);
        end

        fvalues=fval;

        %t_delayVp(ii)=x(5);
        if strcmp(modelType,'full')==1 || strcmp(modelType,'fullTrunc')==1
            if fixedDelay==99
                spillover(ii)=xx(ii,5);
            else
                spillover(ii)=xx(ii,4);
            end
        else
            spillover(ii)=0;
        end
    end


else

    for ii=1:numAzimuthalRegions*numRadialRegions
        kwi(ii) = xx(ii,1);
        kwo(ii)= xx(ii,2);
        %PartCoeff(ii)= xx(ii,2);% this was changed to to replace kwo above (NP 10/10/06)
        est_fb(ii)=fixedVp;
        if fixedDelay==99
            t_delay(ii) = xx(ii,3);
        else
            t_delay(ii)=fixedDelay;   % came from estGlobalDelay
            if fixedVp==99
                est_fb(ii)= xx(ii,3);
            end
        end
        if fixedVp==99 & fixedDelay==99
            est_fb(ii)= xx(ii,4);
        end

        fvalues=fval;

        %t_delayVp(ii)=x(5);
        if strcmp(modelType,'full')==1 || strcmp(modelType,'fullTrunc')==1
            if fixedDelay==99
                spillover(ii)=xx(ii,5);
            else
                spillover(ii)=xx(ii,4);
            end
        else
            spillover(ii)=0;
        end

    end
end

%save tmpParams

if fixedDelay~=99  % ugly hack so won't have param as part of filename
  fixedDelay=0;
end 
% write out in ascii

if exist('seriesNumAIF')
    outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF,flagTimeStamps);
else
    outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp);
end

%outfilename=strcat(outfilename,'.funcfit.txt')
outfilename=strcat(outfilename,'.',modelType,'.txt');

fid=fopen(outfilename,'w');
index=1:length(kwi);
%fprintf(fid,'Created %s ,  delta_t=%f   Order below is kwi, kwo, fv, t_delay, fval\n',date, delta_t);
fprintf(fid,'Created %s ,  delta_t=%f\n',date, delta_t);
%out=[index; kwi./0.5];  % This assumes Flow=K1/E, where E=0.5
% fprintf(fid,'flow  %d    %6.3f\n',out);  % assume E=0.5, kwi=E*F
% fprintf(fid,'\n');

out=[index; kwi];  % EVRD 3/09
fprintf(fid,'Ktrans  %d    %6.3f\n',out);  % assume E=0.5, kwi=E*F
fprintf(fid,'\n');

%fprintf(fid,'flowMean and Std  %6.2f %6.2f , coeff var= %6.2f\n',mean(kwi./0.5), std(kwi./0.5), std(kwi./0.5)/mean(kwi./0.5));
 % above assume E=0.5, kwi=E*F
 fprintf(fid,'flowMean and Std  %6.2f %6.2f , coeff var= %6.2f\n',mean(kwi), std(kwi), std(kwi)/mean(kwi));
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This was added to write out a .mat file of the 1-region (or xx-region etc.) Flow values for
% the select model only.  NP 10/23/06
if strcmp(modelType,'full')==1
    if (numAzimuthalRegions==1)  % this can be changed to any selected numAzimuthalRegions
        if (fixedDelay==0)
            tempFlow=out(2,:);
            FlowFileName=strcat(outpath,'tempFlow.study',int2str(studyNum),'.slice',int2str(sliceNum),'.1region.full_fixedDelay0_fixedVp99_deltaSI.mat');
            save(FlowFileName,'tempFlow');
            %     elseif (fixedDelay==99)   %%%% for noBlood1 model, fixedDelay is always equal to 1.
            %         tempVe=out(2,:);
            %         VeFileName=strcat(outpath,'tempVe.study',int2str(studyNum),'.slice',int2str(sliceNum),'.64region.noBlood1_fixedDelay99_fixedVp99.mat');
            %         save(VeFileName,'tempVe');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hct=0.45; 
out=[index; (1-Hct)*kwi./kwo];   % becuase kwi for arterial
%out=[index; (1-Hct)*PartCoeff]; % this was added to replace kwi/kwo above (NP 10/10/06)
fprintf(fid,'ve    %d     %6.3f\n',out);
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This was added to write out a .mat file of the 8 (or 64 etc.)-region Ve values for
% the noBlood model only.  NP 10/23/06
if strcmp(modelType,'noBlood1')==1
    if (numAzimuthalRegions==8)  % this can be changed to any selected numAzimuthalRegions
        if (fixedDelay==0)
            tempVe=out(2,:);
            VeFileName=strcat(outpath,'tempVe.study',int2str(studyNum),'.slice',int2str(sliceNum),'.8region.noBlood1_fixedDelay0_fixedVp99.mat');
            save(VeFileName,'tempVe');
            %     elseif (fixedDelay==99)   %%%% for noBlood1 model, fixedDelay is always equal to 1.
            %         tempVe=out(2,:);
            %         VeFileName=strcat(outpath,'tempVe.study',int2str(studyNum),'.slice',int2str(sliceNum),'.64region.noBlood1_fixedDelay99_fixedVp99.mat');
            %         save(VeFileName,'tempVe');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out1=(1-Hct)*kwi./kwo;
%out1=(1-Hct)*PartCoeff;  % this was added to replace kwi/kwo above (NP 10/10/06)
%fprintf(fid2,'%6.3f ',out1);





out=[index; est_fb]; 
fprintf(fid,'fv    %d     %6.3f\n',out);
fprintf(fid,'\n');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This was added to write out a .mat file of the 8 (or 64 etc.)-region Vp (est. fv) values for
% % the noBlood model only.  NP 10/23/06
% if strcmp(modelType,'noBlood1')==1
%     if (numAzimuthalRegions==8)  % this can be changed to any selected numAzimuthalRegions
%         if (fixedDelay==0)
%             tempVp=out(2,:);
%             VpFileName=strcat(outpath,'tempVp.study',int2str(studyNum),'.slice',int2str(sliceNum),'.8region.noBlood1_fixedDelay0_fixedVp99.mat');
%             save(VpFileName,'tempVp');
%             %     elseif (fixedDelay==99)   %%%% for noBlood1 model, fixedDelay is always equal to 1.
%             %         tempVe=out(2,:);
%             %         VeFileName=strcat(outpath,'tempVe.study',int2str(studyNum),'.slice',int2str(sliceNum),'.64region.noBlood1_fixedDelay99_fixedVp99.mat');
%             %         save(VeFileName,'tempVe');
%         end
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out1=est_fb;
%fprintf(fid2,'%6.3f ',out1);




out=[index; t_delay]; 
fprintf(fid,'t0    %d     %6.3f\n',out);
fprintf(fid,'\n');

%out=[index; t_delayVp]; 
%fprintf(fid,'t_delayVp    %d     %6.3f\n',out);
%fprintf(fid,'\n');

% out=[index; spillover]; 
% fprintf(fid,'spillover    %d     %6.3f\n',out);
% fprintf(fid,'\n');

%out1=t_delay;
%fprintf(fid2,'%6.3f ',out1);
%fclose(fid2);

out=[index; fvalues]; 
fprintf(fid,'fval  %d     %6.3f\n',out);

fclose(fid);

return

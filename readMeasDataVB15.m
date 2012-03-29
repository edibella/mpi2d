%function kSpace = readMeasDataVB15(filename)
% Program to read measurement data from Siemens MRI scanners
% with IDEA VB15 using Matlab 7 (single value)
%
% function [kSpace] = readMeasDataVB15(filename)
% If input argument 'filename' is not defined
% use GUI to select file with measurements (*.dat)
% kSpace - measured k-space data
%
% Version 1.0
% AUTHOR of MATLAB implementation: Eugene G. Kholmovski, PhD
% UCAIR, Department of Radiology, University of Utah
% ekhoumov@ucair.med.utah.edu


%% GUI to define filename for measurement file
%if nargin == 0
  filename = '';
  if length(filename) == 0
     [temp path] = uigetfile('*.dat','Select File to Read');
     filename = [path temp(1:length(temp)-4)];
     nargin = 0;
  end
%end


%% Definition of OFF and ON constants
OFF = 0;
ON = 1;


%% resetFFTscale should be equal to 1 to reset FFTscale and DataCorrection
%% for each coil to 1
resetFFTscale = OFF;


%% Flags to Read Meas Data from Individual Coil
readOneCoil = OFF;

% coilIndex defines which coil element data will be read
% when readOneCoil = ON
coilIndex = 1;


%% Flags to remove oversampling (OS) in x-direction
%% One of these flags should be equal to 1 to remove OS.
%% removeOS=1 is more efficient as it processes each readout line
%% independently reducing the required memory space to keep all measured
%% data.
removeOS = OFF;  % EVRD 8/27/08, just to try with radial and compare. Not good!  Back to OFF
removeOSafter = OFF;   % note this works in image space, cutting FOV. Not likely good idea for radial


%% transformToImageSpace should be equal to 1 to get image space
%% representation. Take into account that no correction for partial Fourier
%% or parallel imaging k-space undersampling is done.
%% The given version of code only uses FFT operation.
transformToImageSpace = OFF;


%% Flags
readPhaseCorInfo = OFF;
readNavigator = OFF;
readTimeStamp = ON;


%% writeToFile should be equal to 1 to save k-space or image space volume in Matlab file (*.mat)
writeToFile = ON;
if writeToFile == 1
  if transformToImageSpace == 0
     filenameOut = sprintf('%s_Kspace',filename);
  else
     filenameOut = sprintf('%s_imageSpace',filename);
  end
end


%% Useful Parameters
globalHeader = 32;
localHeader = 128;


%% Read Protocol from filename.dat file
if nargin == 1
  fid = fopen(filename,'r','ieee-le');
else
  fid = fopen([filename '.dat'],'r','ieee-le');
end
dataField = fread(fid,1,'int32');
info = fread(fid,dataField-4,'uchar');
info = char(info');
fclose(fid);
longProtocol = info;

textStart = 'MeasYaps';
textEnd = 'Phoenix';
indexStart = strfind(info,textStart) + length(textStart) + 5;
indexEnd = strfind(info,textEnd) - 3;
shortProtocol = info(indexStart:indexEnd);
clear info;
info = shortProtocol;


%% Some Acquisition Parameters from filename.asc
text = 'ucDimension                      = 0x';
t = strfind(info,text) + length(text);
ScanDimension = sscanf(info(t:t+10),'%d');
if ScanDimension == 4
  flag3D = ON;
else
  flag3D = OFF;
end
text = 'sKSpace.lBaseResolution                  = ';
t = strfind(info,text) + length(text);
NxAll = sscanf(info(t:t+10),'%d');
text = 'sKSpace.lPhaseEncodingLines              = ';
t = strfind(info,text) + length(text);
NyAll = sscanf(info(t:t+10),'%d');
text = 'sKSpace.dPhaseOversamplingForDialog      = ';
t = strfind(info,text) + length(text);
if isempty(t) == 1
   OSfactorPE = 1;
else
   OSfactorPE = 1 + sscanf(info(t:t+10),'%f');
end
text = 'sKSpace.lPartitions                      = ';
t = strfind(info,text) + length(text);
NzAll = sscanf(info(t:t+10),'%d');
text = 'sKSpace.dSliceOversamplingForDialog      = ';
t = strfind(info,text) + length(text);
if isempty(t) == 1
   OSfactor3D = 1;
else
   OSfactor3D = 1 + sscanf(info(t:t+10),'%f');
end
text = 'sKSpace.dPhaseResolution                 = ';
t = strfind(info,text) + length(text);
phaseResolution = sscanf(info(t:t+10),'%f');
text = 'sKSpace.dSliceResolution                 = ';
t = strfind(info,text) + length(text);
sliceResolution = sscanf(info(t:t+10),'%f');
text = 'sSliceArray.lSize                        = ';
t = strfind(info,text) + length(text);
Nsl = sscanf(info(t:t+10),'%d');
text = 'sSliceArray.lConc                        = ';
t = strfind(info,text) + length(text);
Nconc = sscanf(info(t:t+10),'%d');
text = '.lRxChannelConnected = ';
t = strfind(info,text) + length(text);
Nc = length(t);
text = 'AdjustSeq%/AdjCoilSensSeq';
t = strfind(info,text) + length(text);
if isempty(t) == 0
  Nc = Nc-1;
end
text = 'lContrasts                               = ';
t = strfind(info,text) + length(text);
if isempty(t) == 1
  nContrast = 1;
else
  nContrast = sscanf(info(t:t+10),'%d');
end
text = 'lSets                                    = ';
t = strfind(info,text) + length(text);
if isempty(t) == 1
  nSet = 1;
else
  nSet = sscanf(info(t:t+10),'%d');
end
text = 'lAverages                                = ';
t = strfind(info,text) + length(text);
if isempty(t) == 1
  nAverage = 1;
else
  nAverage = sscanf(info(t:t+10),'%d');
end
text = 'lRepetitions                             = ';
t = strfind(info,text) + length(text);
if isempty(t) == 1
   nRepetition = 1;
else
   nRepetition = sscanf(info(t:t+10),'%d') + 1;
end
text = 'sPhysioImaging.lPhases                   = ';
t = strfind(info,text) + length(text);
if isempty(t) == 1
   nPhase = 1;
else
   nPhase = sscanf(info(t:t+10),'%d');
end
text = 'sKSpace.PhasePartialFourier            = 0x';
t = strfind(info,text) + length(text);
if isempty(t) == 1
  fractionFlag = 10;
else
  fractionFlag = str2num(info(t:t+1));
end
switch(fractionFlag)
   case(10)
       fractionY = 1.0;
   case(8)
       fractionY = 7/8;
   case(4)
       fractionY = 0.75;
   case(2)
       fractionY = 5/8;
   case(1)
       fractionY = 0.5;
end
text = 'sKSpace.SlicePartialFourier            = 0x';
t = strfind(info,text) + length(text);
if isempty(t) == 1
  fractionFlag = 10;
else
  fractionFlag = str2num(info(t:t+1));
end
switch(fractionFlag)
   case(10)
       fractionZ = 1.0;
   case(8)
       fractionZ = 7/8;
   case(4)
       fractionZ = 0.75;
   case(2)
       fractionZ = 5/8;
   case(1)
       fractionZ = 0.5;
end
text = 'sKSpace.dSeqPhasePartialFourierForSNR    = ';
t = strfind(info,text) + length(text);
if isempty(t) == 1
  phasePartialFourierForSNR = 1.0;
else
  phasePartialFourierForSNR = sscanf(info(t:t+10),'%f');
end

text = 'sFastImaging.lEPIFactor                  = ';
t = strfind(info,text) + length(text);
EPIFactor = sscanf(info(t:t+10),'%d');
text = 'sFastImaging.lTurboFactor                = ';
t = strfind(info,text) + length(text);
if isempty(t) == 1
   turboFactor = 1;
else
   turboFactor = sscanf(info(t:t+10),'%d');
end


%% iPAT parameters from filename.asc
text = 'sPat.ucPATMode                           = 0x';
t = strfind(info,text) + length(text);
PATMode = str2num(info(t:t+1));
text = 'sPat.ucRefScanMode                       = 0x';
t = strfind(info,text) + length(text);
PATRefScanMode = str2num(info(t:t+1));
text = 'sPat.lAccelFactPE                        = ';
t = strfind(info,text) + length(text);
AccelFactorPE = sscanf(info(t:t+10),'%d');
text = 'sPat.lAccelFact3D                        = ';
t = strfind(info,text) + length(text);
AccelFactor3D = sscanf(info(t:t+10),'%d');
if AccelFactorPE == 1
   nRefLinesPE = 0;
else
   text = 'sPat.lRefLinesPE                         = ';
   t = strfind(info,text) + length(text);
   nRefLinesPE = sscanf(info(t:t+10),'%d');
end
if AccelFactor3D == 1
   nRefLines3D = 0;
else
   text = 'sPat.lRefLines3D                         = ';
   t = strfind(info,text) + length(text);
   nRefLines3D = sscanf(info(t:t+10),'%d');
end


clear info;
info = longProtocol;
text = 'ParamLong."iMaxNoOfRxChannels';
t = strfind(info,text) + length(text);
if isempty(t) == 0
  ss = sscanf(info(t:t+35),'%s');
  Nct = str2num(ss(isstrprop(ss,'digit')));
  if Nct < Nc
     Nc = Nct;
  end
end

text = 'ParamLong."lNoOfPhaseCorrScans';
t = strfind(info,text) + length(text);
if isempty(t) == 0
  ss = sscanf(info(t:t+35),'%s');
  if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
     nPhCorrScan = str2num(ss(isstrprop(ss,'digit')));
  end
end
if turboFactor > 1
  nPhCorEcho = 1;
  nPhCorScan = 1;
end
if EPIFactor > 1
  nPhCorScan = 1;
  nPhCorEcho = 3;
end

if AccelFactorPE == 1
   FirstFourierLine = 1;
   FirstRefLine = 1;
else
   text = 'ParamLong."lFirstFourierLine';
   t = strfind(info,text) + length(text);
   FirstFourierLine = 0;
   if isempty(t) == 0
      ss = sscanf(info(t:t+35),'%s');
      if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
         FirstFourierLine = str2num(ss(isstrprop(ss,'digit')));
      end
   end
   FirstFourierLine = FirstFourierLine + 1;
   text = 'ParamLong."lFirstRefLine';
   t = strfind(info,text) + length(text);
   FirstRefLine = 0;
   if isempty(t) == 0
      ss = sscanf(info(t:t+35),'%s');
      if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
         FirstRefLine = str2num(ss(isstrprop(ss,'digit')));
      end
   end
   FirstRefLine = FirstRefLine + 1;
end
if AccelFactor3D == 1
   FirstFourierPar = 1;
   FirstRefPartition = 1;
else
   text = 'ParamLong."lFirstFourierPartition';
   t = strfind(info,text) + length(text);
   FirstFourierPartition = 0;
   if isempty(t) == 0
      ss = sscanf(info(t:t+35),'%s');
      if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
         FirstFourierPartition = str2num(ss(isstrprop(ss,'digit')));
      end
   end
   FirstFourierPartition = FirstFourierPartition + 1;
   text = 'ParamLong."lFirstRefPartition';
   FirstRefPartition = 0;
   t = strfind(info,text) + length(text);
   if isempty(t) == 0
      ss = sscanf(info(t:t+35),'%s');
      if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
         FirstRefPartition = str2num(ss(isstrprop(ss,'digit')));
         if FirstRefPartition < 0
            FirstRefPartition = 0;
         end
      end
   end
   FirstRefPartition = FirstRefPartition + 1;
end

%% YAPS Parameters from filename.asc
text = 'iNoOfFourierColumns';
t = strfind(info,text) + length(text);
ss = sscanf(info(t:t+35),'%s');
NxOS = str2num(ss(isstrprop(ss,'digit')));

text = 'flReadoutOSFactor';
t = strfind(info,text) + length(text);
OSfactorRO = sscanf(info(t:t+10),'%f');
OSfactorRO = 2;
Nx = round(NxOS/OSfactorRO);

text = 'NoOfFourierLines';
t = strfind(info,text) + length(text);
ss = sscanf(info(t:t+35),'%s');
Ny = str2num(ss(isstrprop(ss,'digit')));
text = 'iNoOfFourierPartitions';
t = strfind(info,text) + length(text);
ss = sscanf(info(t:t+35),'%s');
Nz = str2num(ss(isstrprop(ss,'digit')));
text = 'iRoFTLength';
t = strfind(info,text) + length(text);
ss = sscanf(info(t:t+35),'%s');
NxRecon = str2num(ss(isstrprop(ss,'digit')));
text = 'iPEFTLength';
t = strfind(info,text) + length(text);
ss = sscanf(info(t:t+35),'%s');
NyRecon = str2num(ss(isstrprop(ss,'digit')));
text = 'i3DFTLength';
t = strfind(info,text) + length(text);
ss = sscanf(info(t:t+35),'%s');
NzRecon = str2num(ss(isstrprop(ss,'digit')));
if Nz == 1
  NzRecon = 1;
end
text = 'ParamLong."ushSlicePerConcat';
t = strfind(info,text) + length(text);
ss = sscanf(info(t:t+35),'%s');
NslicePerConcat = str2num(ss(isstrprop(ss,'digit')));


%% Partial Fourier Mode and Parameters from filename.asc
text = 'lPCAlgorithm';
t = strfind(info,text) + length(text);
ss = sscanf(info(t:t+35),'%s');
PCAlgorithm = str2num(ss(isstrprop(ss,'digit')));
text = 'ParamLong."lNoOfPhaseCorrColumns';
t = strfind(info,text) + length(text);
NoOfPhaseCorrColumns = Nx/2;
if isempty(t) == 0
  ss = sscanf(info(t:t+35),'%s');
  if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
     NoOfPhaseCorrColumns = str2num(ss(isstrprop(ss,'digit')));
  end
end
text = 'ParamLong."lNoOfPhaseCorrLines';
t = strfind(info,text) + length(text);
NoOfPhaseCorrLines = NyAll/2;
if isempty(t) == 0
  ss = sscanf(info(t:t+35),'%s');
  if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
     NoOfPhaseCorrLines = str2num(ss(isstrprop(ss,'digit')));
  end
end
text = 'ParamLong."lNoOfPhaseCorrPartitions';
t = strfind(info,text) + length(text);
NoOfPhaseCorrPartitions = NzAll/2;
if isempty(t) == 0
  ss = sscanf(info(t:t+35),'%s');
  if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
     NoOfPhaseCorrPartitions = str2num(ss(isstrprop(ss,'digit')));
  end
end
text = 'ParamLong."lColSlopeLength';
t = strfind(info,text) + length(text);
ColSlopeLength = Nx/4;
if isempty(t) == 0
  ss = sscanf(info(t:t+35),'%s');
  if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
     ColSlopeLength = str2num(ss(isstrprop(ss,'digit')));
  end
end
text = 'ParamLong."lLinSlopeLength';
t = strfind(info,text) + length(text);
ColSlopeLength = NyAll/4;
if isempty(t) == 0
  ss = sscanf(info(t:t+35),'%s');
  if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
     LinSlopeLength = str2num(ss(isstrprop(ss,'digit')));
  end
end
text = 'ParamLong."lParSlopeLength';
t = strfind(info,text) + length(text);
ParSlopeLength = NzAll/4;
if isempty(t) == 0
  ss = sscanf(info(t:t+35),'%s');
  if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
     ParSlopeLength = str2num(ss(isstrprop(ss,'digit')));
  end
end


clear info;
info = shortProtocol;
%% Raw Data Correction Factors
CorrFactor = ones(Nc,1);
if ON == 0
for c=1:Nc
  text = 'axRawDataCorrectionFactor[0][';
  text = sprintf('%s%d',text,c-1);
  if c < 11
     text = [text,'].dRe      ='];
  else
     text = [text,'].dRe     ='];
  end
  t = strfind(info,text) + length(text);
  if isempty(t) == 1
     CorrFactor(c) = 1;
  else
     CorrFactor(c) = sscanf(info(t:t+20),'%f');
  end
  text = 'axRawDataCorrectionFactor[0][';
  text = sprintf('%s%d',text,c-1);

  if c < 11
     text = [text,'].dIm      ='];
  else
     text = [text,'].dIm     ='];
  end
  t = strfind(info,text) + length(text);
  if isempty(t) == 0
     CorrFactor(c) = CorrFactor(c) + i*sscanf(info(t:t+20),'%f');
  end
end
end


%% FFT Correction Factors
FFTCorrFactor = ones(Nc,1);
for c=1:Nc
  text = 'aFFT_SCALE[';
  text = sprintf('%s%d',text,c-1);
  if c < 11
    text = [text,'].flFactor = '];
  else
    text = [text,'].flFactor = '];
  end
  
 % text = [text,'].flFactor'];   % evrd 4/08  to read brain 12 coil
  t = strfind(info,text) + length(text);
 % t=t+3;   % evrd to get back the space = space
 % keyboard
  FFTCorrFactor(c) = sscanf(info(t:t+20),'%f');
 % FFTCorrFactor(c) = sscanf(info(t+10:t+18),'%f');
  
end
if resetFFTscale == 1
  FFTCorrFactor = ones(Nc,1);
end


%% For PC Angio
Nset = 1;
text = 'sAngio.ucPCFlowMode                      = 0x';
t = strfind(info,text) + length(text);
if isempty(t) == 0
  PCMRAFlag = sscanf(info(t:t+1),'%d');
else
  PCMRAFlag = 0;
end
if PCMRAFlag == 1
  text = 'sAngio.sFlowArray.lSize                  = ';
  t = strfind(info,text) + length(text);
  if isempty(t) == 0
     Nset = sscanf(info(t:t+2),'%d') + 1;
  end
end

clear info;


%% Recalculation of partial Fourier factors and EPI/turbo factor for EPI and TSE
fractionPE = Ny/NyAll;
fraction3D = Nz/NzAll;
EPIFactor = round(fractionPE*EPIFactor);
nEPITrain = Ny/EPIFactor;
turboFactor = round(fractionPE*turboFactor);
nTSETrain = Ny/turboFactor;


Nc0 = Nc;
if readOneCoil == 1
  Nc = 1;
else
  Nc = Nc0;
end

%% Calculation of the number of valid k-space readouts and k-space data matrix dimensions
if PATMode == 1
   nReadout = nAverage*nPhase*nRepetition*nContrast*Nsl*Nz*Nc*Ny;
elseif PATMode == 2 & PATRefScanMode == 2
   if mod(Ny,2) == 1
      NyPAT = (Ny-1+nRefLinesPE*(AccelFactorPE-1))/AccelFactorPE;
   else
      NyPAT = floor((Ny+nRefLinesPE*(AccelFactorPE-1))/AccelFactorPE);
   end
   nReadout = nAverage*nPhase*nRepetition*nContrast*Nsl*Nz*Nc*NyPAT;
end

if removeOS == 1
   kSpace = zeros(nAverage, nPhase, nRepetition, nContrast, Nsl, Nc, Nz, Nx, Ny, 'single');
else
   kSpace = zeros(nAverage, nPhase, nRepetition, nContrast, Nsl, Nc, Nz, NxOS, Ny, 'single');
end
if readPhaseCorInfo == 1 & nPhCorScan > 0
   kPhaseCor = zeros(nPhCorScan, nPhCorEcho, Nsl, nRepetition, Nc, NxOS, 'single');
end
if readNavigator == 1
   nNavigator = nAverage*nPhase*nRepetition*nContrast*Nsl*Nz*Nc*nEPITrain*nNavEK;
   kNavigator = zeros(nAverage, nPhase, nRepetition, nContrast*nNavEK, Nsl, Nc, Nz, nEPITrain, NxOS, 'single');
end
if readTimeStamp == 1
  nTimeStamp = nAverage*nPhase*nRepetition*nContrast*Nz;
  timeStamp = zeros(nAverage, nPhase, nRepetition, nContrast*nNavEK, Nz, 'single');
end

%keyboard;

%%% DATA READOUT & REORDERING
%% Read k-space data from filename.dat
noiseMeasCounter = 0;
noiseMeas = zeros(NxOS,Nc);
navigatorPrep = 0;
LineN = -1;
Ndr = 1;
for r=1:Ndr
  xCoil = r-1
  if nargin == 1
     fid = fopen(filename,'r','ieee-le');
  else
     fid = fopen([filename '.dat'],'r','ieee-le');
  end
  readFlag = 1;
  skipField = 0;
  count = 0;
  countNavigator = 0;
  navigatorDataON = 0;
  temp1 = zeros(NxOS,1);
  dataField = fread(fid,1,'int32');
  st = fseek(fid,dataField,-1);
  while (readFlag == 1)
      if readTimeStamp  == 1
         st = fseek(fid,12,0);
         timeS = fread(fid,1,'uint32');
         st = fseek(fid,4,0);
      else
         st = fseek(fid,20,0);
      end
      evalMask1 = fread(fid,1,'uint32');
      evalMask2 = fread(fid,1,'uint32');
      flag = 33-find(dec2base(evalMask1,2,32)  == '1');
      Nxr = fread(fid,1,'uint16');
      Ncr = fread(fid,1,'uint16');
      Line = fread(fid,1,'uint16');
      Acquisition = fread(fid,1,'uint16');
      Slice = fread(fid,1,'uint16');
      Partition = fread(fid,1,'uint16');
      Echo = fread(fid,1,'uint16');
      Phase = fread(fid,1,'uint16');
      Repetition = fread(fid,1,'uint16');
      Set = fread(fid,1,'uint16');
      st = fseek(fid,12,0);
      CutOffDataPre = fread(fid,1,'uint16');
      CutOffDataPost = fread(fid,1,'uint16');
      KSpaceCentreColumn = fread(fid,1,'uint16');
      CoilMode = fread(fid,1,'uint16');
      ReadOutOffCentre = fread(fid,1,'float32');
      st = fseek(fid,4,0);
      KSpaceCentreLineNo = fread(fid,1,'uint16');
      KSpaceCentrePartitionNo = fread(fid,1,'uint16');
      st = fseek(fid,44,0);
      Channel = fread(fid,1,'uint16');
      st = fseek(fid,2,0);

      if isempty(find(flag == 1)) == 0
         break;
      end


      if isempty(find(flag == 2)) == 0 | isempty(find(flag == 22)) == 0 | isempty(find(flag == 26)) == 0
         if isempty(find(flag == 22)) == 0 & readPhaseCorInfo == 0
            skipField = nPhCorScan*nPhCorEcho*Ncr*(localHeader+8*Nxr)-localHeader;
         end
         if isempty(find(flag == 22)) == 0 & readPhaseCorInfo == 1
            skipField = -localHeader;
            st = fseek(fid,skipField,0);
            skipField = 0;
            for m = 1:nPhCorScan*nPhCorEcho*Ncr
               infoMDH_TimeStamp = readMDH_TimeStamp_VB13(fid);
               temp = fread(fid,2*Nxr,'float32');
               if isempty(find(flag == 25)) == 0
                  temp(1:2:end) = flipud(temp(1:2:end));
                  temp(2:2:end) = flipud(temp(2:2:end));
               end
               if CutOffDataPre > 0
                  temp(1:2*CutOffDataPre) = 0;
               end
               if CutOffDataPost > 0
                  temp(end-2*CutOffDataPost+1:end) = 0;
               end
               kPhaseCor(Echo+1, ceil(m/Ncr), Slice+1, Repetition+1, Channel+1, :) = single(temp(1:2:end)+i*temp(2:2:end));
            end
          end

          if isempty(find(flag == 2)) == 0 & readNavigator == 0
             skipField = Ncr*(localHeader+8*Nxr)-localHeader;
          end
          if isempty(find(flag == 2)) == 0 & readNavigator == 1
                if countNavigator == 0 & navigatorPrep == 0
                   kNavigator = single(zeros(Nxr, Ncr, nContrast*nNavEK, nEPITrain, Nz, Nsl, nAverage, nPhase, nRepetition));
                   kNavigatorTemp = zeros(Nxr, Ncr, nContrast*nNavEK);
                   navigatorPrep = 1;
                end
                skipField = -localHeader;
                st = fseek(fid,skipField,0);
                skipField = 0;
                for m=1:nNavEK*Ncr
                   infoMDH_TimeStamp = readMDH_TimeStamp_VB13(fid);
                   temp = fread(fid,2*Nxr,'float32');
                   if isempty(find(flag == 25)) == 0
                      temp(1:2:end) = flipud(temp(1:2:end));
                      temp(2:2:end) = flipud(temp(2:2:end));
                   end
                   if CutOffDataPre > 0
                      temp(1:2*CutOffDataPre) = 0;
                   end
                   if CutOffDataPost > 0
                      temp(end-2*CutOffDataPost+1:end) = 0;
                   end
                   kNavigatorTemp(:,Channel+1,Echo+1) = single(temp(1:2:end)+i*temp(2:2:end));
                end
                navigatorDataON = 1;
          end

          if isempty(find(flag == 26)) == 0
             temp = fread(fid,2*Nxr,'float32');
             if isempty(find(flag == 25)) == 0
                temp(1:2:end) = flipud(temp(1:2:end));
                temp(2:2:end) = flipud(temp(2:2:end));
             end
             noiseMeas(:,Channel+1) = temp(1:2:end) + i*temp(2:2:end);
             skipField = 0;
          end
          st = fseek(fid,skipField,0);

      else
          temp = fread(fid,2*Nxr,'float32');
          if isempty(find(flag == 25)) == 0
             temp(1:2:end) = flipud(temp(1:2:end));
             temp(2:2:end) = flipud(temp(2:2:end));
          end
          if CutOffDataPre > 0
             temp(1:2*CutOffDataPre) = 0;
          end
          if CutOffDataPost > 0
             temp(end-2*CutOffDataPost+1:end) = 0;
          end
          if isempty(find(flag == 11)) == 0
              temp = CorrFactor(Channel+1)*temp;
          end
          temp = FFTCorrFactor(Channel+1)*temp;


          if readOneCoil == 0
             if removeOS == 1
                temp1(end-Nxr+1:end) = temp(1:2:end)+i*temp(2:2:end);
                tempX = fftshift(fft(fftshift(temp1)));
                tempK = fftshift(ifft(fftshift(tempX(round((NxOS-Nx)/2)+1:Nx+round((NxOS-Nx)/2)))));
                kSpace(Acquisition+1, Phase+1, Repetition+1, Echo+1, Slice+1, Channel+1, Partition+1, :, Line+1) = single(tempK);
             else
                kSpace(Acquisition+1, Phase+1, Repetition+1, Echo+1, Slice+1, Channel+1, Partition+1, end-Nxr+1:end, Line+1) = single(temp(1:2:end)+i*temp(2:2:end));
             end
             count = count+1;
             if mod(count,10) == 0
                 count;
             end
          elseif readOneCoil == 1 & Channel+1 == coilIndex
             if removeOS == 1
                temp1(end-Nxr+1:end) = temp(1:2:end)+i*temp(2:2:end);
                tempX = fftshift(fft(fftshift(temp1)));
                tempK = fftshift(ifft(fftshift(tempX(round((NxOS-Nx)/2)+1:Nx+round((NxOS-Nx)/2)))));
                kSpace(Acquisition+1, Phase+1, Repetition+1, Echo+1, Slice+1, 1, Partition+1, :, Line+1) = single(tempK);
             else
                kSpace(Acquisition+1, Phase+1, Repetition+1, Echo+1, Slice+1, 1, Partition+1, end-Nxr+1:end, Line+1) = single(temp(1:2:end)+i*temp(2:2:end));
             end
             count = count+1;
          end

          if readTimeStamp == 1 & Channel == 0 & navigatorDataON == 1
             EPITrain = mod(countNavigator, nEPITrain);
             timeStamp(Echo+1, EPITrain+1, Partition+1, Slice+1, Acquisition+1, Phase+1, Repetition+1) = single(0.0025*timeS);
          end

          if readNavigator == 1 & Channel == 0 & navigatorDataON == 1
             EPITrain = mod(countNavigator, nEPITrain);
             kNavigator(:, :, :, EPITrain+1, Partition+1, Slice+1, Acquisition+1, Phase+1, Repetition+1) = single(kNavigatorTemp);
             navigatorDataON = 0;
             countNavigator = countNavigator+1;
          end
      end
      if isempty(find(flag == 1)) == 0
         break;
      end
  end
   fclose(fid);
   clear temp temp1 tempX tempK;
end
kSpace = squeeze(kSpace);


if ndims(kSpace) == 3
  kSpace = permute(kSpace,[2 3 1]);
end
if ndims(kSpace) == 4
  if flag3D == 0
     kSpace = permute(kSpace,[3 4 1 2]);
  else
     kSpace = permute(kSpace,[3 4 2 1]);
  end
end
if ndims(kSpace) == 5
  kSpace = permute(kSpace,[4 5 3 1 2]);
end


DataName = 'kSpace';
if transformToImageSpace == 1
  imSpace = single(NxOS*fftshift(ifft(fftshift(kSpace,1),[],1),1));
  imSpace = Ny*fftshift(ifft(fftshift(imSpace,2),[],2),2);
  if flag3D == 1
     imSpace = Nz*fftshift(ifft(fftshift(imSpace,3),[],3),3);
  end
  if removeOSafter == 1 && removeOS == 0
     if ndims(imSpace) == 2
        imSpace(1:NxOS/4,:) = [];
        imSpace(end-NxOS/4+1:end,:) = [];
     end
     if ndims(imSpace) == 3
        imSpace(1:NxOS/4,:,:) = [];
        imSpace(end-NxOS/4+1:end,:,:) = [];
     end
     if ndims(imSpace) == 4
        imSpace(1:NxOS/4,:,:,:) = [];
        imSpace(end-NxOS/4+1:end,:,:,:) = [];
     end
     if ndims(imSpace) == 5
        imSpace(1:NxOS/4,:,:,:,:) = [];
        imSpace(end-NxOS/4+1:end,:,:,:,:) = [];
     end
  end
  DataName = 'imSpace';
end

if writeToFile == 1
  save(filenameOut, DataName,'timeStamp');
end

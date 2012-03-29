PatientName = 'P030909';

VerTrio = '';
cd('/v/raid1/npack/MRIdata/Cardiac');
if(exist(['Verio/' PatientName]))
    VerTrio = 'Verio';
elseif(exist(['Trio/' PatientName]))
    VerTrio = 'Trio';
else
    disp(['Could not find ' PatientName ' in Nate''s MRI folder']);
    return;
end

cd([VerTrio '/' PatientName]);
NatesPatientFolder = pwd;
cd(['/v/raid1/bmatthew/MRIdata/Cardiac/' VerTrio]);
if(exist(PatientName))
    disp(['Brian already has ' PatientName]);
    return;
else
    mkdir(PatientName);
    cd(PatientName);
    PatientFolder = pwd;
    mkdir('Processing'); cd('Processing'); processingFolder = pwd; cd('..');
    mkdir('DicomData') ; cd('DicomData');  DicomFolder = pwd;      cd('..');
    mkdir('ReconData');
    mkdir('ReconData/mat_files');
end

if(exist('DicomData') == 7)
    cd(PatientFolder);
    DoThisFolder = ['/v/raid1/npack/MRIdata/Cardiac/' VerTrio '/' PatientName];
    command = ['mpi2d stg=0.1 DoThisFolder=''' DoThisFolder ''''];
    eval(command);
end
cd(NatesPatientFolder);
if(exist('RawData') == 7)
    if(exist([PatientFolder '/RawData']) ~= 7)
        disp(['Copying RawData folder from ' NatesPatientFolder]);
        copyfile('RawData',[PatientFolder '/RawData']);
    end
end

if(exist('ReconData') == 7)
    cd('ReconData');
    recons = dir('*');
    for i=3:length(recons)
        
    end
end
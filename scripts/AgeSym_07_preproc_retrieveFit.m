function AgeSym_07_preproc_retrieveFit(cohort)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Purpose: retrieve AgexVertex csv file for signficant 
    %          asymmetry trajectories and hemi effects        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    base='/Users/jamesroe/Dropbox/OpenScienceFramework/AgeSym';
    cohort='CamCan'
    resdir=[base filesep 'results/reproduceClustering/' cohort];
    labelname=[resdir filesep 'vtxfsav5.csv'];
    
    addpath(genpath([getenv('FREESURFER_HOME') filesep 'matlab']));

    % read label
    a = fopen(labelname);
    vtx = textscan(a,'%s');
    lab = str2double(vtx{1}) + 1;


    cd(resdir)
    
    %fit
    filename = 'fitFitfsav5.csv';
    if ~exist(filename)
    
        [hdvol, hdM, hdmr_parms, hdvolsz] = load_mgh('m1004dfit.fsav5.mgh');
    
        D = []
        for i = 1:100
           disp(i);
           D(:,i) = hdvol(lab(:,1),1,i);
        end
        disp(filename)
        csvwrite(filename,D);
    end
    
    
    %Hemisphere effects
    filename = 'mapHCoef.fsav5.csv'; 
    if ~exist(filename)
        
        [vol, M, mr_parms, volsz] = load_mgh('mapHCoef.fsav5.mgh');
    
        D = []
        D = vol(lab(:,1));
        disp(filename)
        csvwrite(filename,D);
    end
    quit()
end


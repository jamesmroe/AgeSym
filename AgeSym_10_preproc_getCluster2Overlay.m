function AgeSym_10_preproc_getCluster2Overlay(cohort)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Purpose: make fsaverage5 surface overlay of clustering results 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    base='/Users/jamesroe/Dropbox/OpenScienceFramework/AgeSym';
    cohort='LCBC'
    resdir=[base filesep 'results/reproduceClustering/' cohort];
    
    
    % read vertices of label
    addpath(genpath([getenv('FREESURFER_HOME') filesep 'matlab']));
    labelname = 'vtxfsav5.csv';
    a = fopen(labelname);
    vtx = textscan(a,'%s');
    lab = str2double(vtx{1}) + 1; %ensure 1-based index
    
	
    % set resultsdir
	plotdir=[resdir filesep 'Results'];
    

	% open image as mgh - template
    cd(resdir)
	mm = dir('finalSigFDR30.LCBC.fsav5.mgh');
	[vol, M, mr_parms, volsz] = load_mgh(mm(1).name);
	vol(:) = 0;


	%get cluster solutions 
	cd(plotdir);
	sols=dir('s.order_*.txt');

	for i = 1:length(sols)
        [~, name,ext] = fileparts(sols(i).name); 
        if ~exist([name '.mgh'])
            disp(name)
            fid = fopen(sols(i).name); 
            A = fscanf(fid, '%f');
            vol(lab(:,1)) = A;
            save_mgh(vol, [plotdir filesep name, '.mgh'], M, mr_parms);
        end
	end 

	quit()
end

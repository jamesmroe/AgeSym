function AgeSym_05_getStats2Overlay(datadir, outdir)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Purpose: make surface overlays from GAMM results 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % read vertices of label
    addpath(genpath([getenv('FREESURFER_HOME') filesep 'matlab']));
    lab = fs_read_label([getenv('SUBJECTS_DIR') '/40tvs8ref_sym_20/label/lh.cortex.label'])';
    
    % open first image as template
	cd(datadir);
	mm = dir('*mgh');
	[vol, M, mr_parms, volsz] = load_mgh(mm(1).name);
	vol(:) = 0;

	% get maps 
	cd(outdir);
    maps=dir('map*.csv');

	for i = 1:length(maps)
        [~, name,ext] = fileparts(maps(i).name); 
        filename = [outdir filesep name, '.mgh'];
        if ~exist(filename)
            disp(name)
            fid = fopen(maps(i).name); 
            A = csvread(maps(i).name);
            vol(lab(:,1)) = A;
            save_mgh(vol, filename, M, mr_parms);
        end
    end
    
    f = 'fitValDiff.csv'
    fid = fopen(f); 
    D = csvread(f);
    
    %load 4D mgh with 100 frames
    [hdvol, hdM, hdmr_parms, hdvolsz] = load_mgh('1004d.mgh');
    hdvol(:,:,:,:) = 0;
 
    %map trajectories to surface
    for i = 1:100
       hdvol(lab(:,1),1,i) = D(:,i);
    end
    save_mgh(hdvol, '1004dfit.mgh', hdM, hdmr_parms);
    
end

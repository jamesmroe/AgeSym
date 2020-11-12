function AgeSym_02_getLabelxSub_Matrix(datadir, outdir, subset, j)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Purpose: write vertex X subject matrix split across N files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    subset=str2double(subset)
    j=str2double(j)
    addpath(genpath([getenv('FREESURFER_HOME') filesep 'matlab']));
    
    
    % read vertices of label
    lab=fs_read_label([getenv('SUBJECTS_DIR') '/40tvs8ref_sym_20/label/lh.cortex.label'])';
    
    splitdir=[outdir filesep 'split.data']
    if ~exist(splitdir)
        mkdir(splitdir);
    end

	% get all raw data
	cd(datadir)
	mm = dir('*mgh');
    subs = length(mm);
    nvtx = length(lab(:,1))
         
    
    % extract Vertex X Sub matrix across N files
    N = ceil(nvtx/subset) 
    C = [];
    disp([num2str(j) '/' num2str(N)])
    if (j == N)
        vtxindex = 1+(j-1)*subset:1:nvtx;
    else
        vtxindex = 1+(j-1)*subset:1:j*subset;
    end
    index = lab(vtxindex);
    disp(['indexing ' num2str(vtxindex(1)) ':' num2str(vtxindex(end))])
    reverseStr = '';
    for i = 1:subs

       % Display progress
       percentDone = 100 * i / subs;
       msg = sprintf('Percent done: %3.1f', percentDone);
       fprintf([reverseStr, msg]);
       reverseStr = repmat(sprintf('\b'), 1, length(msg));

       % Read subject volume
       [vol, M, mr_parms, volsz] = load_mgh(mm(i).name);
       C(:,i) = vol(index);
    end
    disp(size(C))
    %csvwrite([splitdir filesep 'LabelxSub_matrix' num2str(j,'%03.f') '.csv'],C)

end

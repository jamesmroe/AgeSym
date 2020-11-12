%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Compute Dice similarity between  
% clustering results in LCBC and replication samples
% inputs are cortical labels in freesurfer format found in 
% path/to/publishedMaps/labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath([getenv('FREESURFER_HOME') filesep 'matlab']));

clear all
wdir='/Users/jamesroe/Dropbox/OpenScienceFramework/AgeSym/publishedMaps';
cref = 'LCBC'
clabs={'CamCan','BaseII','Betula','DLBS'}
cd(wdir)
for c=1:numel(clabs)
   if (c==1)
       clab='CamCan'
   end
   if (c==2)
       clab='BaseII'
   end
   if (c==3)
       clab='Betula'
   end
   if (c==4)
       clab='DLBS'
   end
   
    for j=1:3 
        
        labref = [wdir '/labels/cl' num2str(j) '-LCBC.label'];
        lab = [wdir '/labels/cl' num2str(j) '-' clab '.label'];
        
        labels = {labref,lab};
        
        C = [];
        C.vector = []; 
        for l=1:numel(labels)
            x = fs_read_label(labels{l})'; %vtxs: Indices of the vertices (1-based). 
            C(l).vector = x; %store vertices 
        end


        % Dice coefficient %
        %%%%%%%%%%%%%%%%%%%%
        %Formula  Dice: S = 2*(A && B) / (A + B);

        if (j==1)
            ss = []; % store dices
            pervtx = [];
        end
        i=1
        k=2
        A = C(i).vector;
        B = C(k).vector;

        S = 2*(numel(intersect(A,B)))/(numel(A) + numel(B));
        [filepath,name1,ext] = fileparts(labels{i})
        [filepath,name2,ext] = fileparts(labels{k})
        fprintf('For labels %s and %s Dice coefficient is: \n%6.4f\n',[cref ' ' name1], [clab ' ' name2], S);

        ss(j) = S;
        
        %percentage of overlapping vertices
        pervtx(j) = numel(intersect(A,B))/numel(B)
                
    end
    
    if (c==1)
        overlaps=[];
        per=[];
    end
    overlaps(c,:)=ss
    perc(c,:) = pervtx
end
%mean across each cluster solution per cohort
%empirical dice values
mm=[mean(overlaps(1,:)),
mean(overlaps(2,:)),
mean(overlaps(3,:)),
mean(overlaps(4,:))]


%%
%permute significance of overlap with LCBC results
volume = [wdir filesep 'Fig3C-fsaverage5-clusteredTrajectories-LCBC.mgh'];
[vol, M, mr_parms, volsz] = load_mgh(volume);
lab=find(vol~=0);
v=vol(lab());


%LCBC reference labels
labref1 = [wdir '/labels/cl1-' cref '.label'];
labref2 = [wdir '/labels/cl2-' cref '.label'];
labref3 = [wdir '/labels/cl3-' cref '.label'];
labels={labref1,labref2,labref3};
C = [];
C.vector = []; 
for l=1:numel(labels)
    x = fs_read_label(labels{l})'; 
    C(l).vector = x(:,1); %store vertices of each solution
end

for c=1:numel(clabs)
    if (c==1)
       clab='CamCan'
    end
    if (c==2)
       clab='BaseII'
    end
    if (c==3)
       clab='Betula'
    end
    if (c==4)
       clab='DLBS'
    end


    %set outputs
    per=10000;
    SSS=[];
    mS=[];
    pS=[];
    reverseStr = '';
    for j=1:per

        % Display the progress
        percentDone = 100 * j / per;
        msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));

        p=v(randperm(length(v)));
        %which vtxs within label
        v1=find(p==1); 
        v2=find(p==2);
        v3=find(p==3);

        %solution 1
        A = C(1).vector; %LCBC vtxs
        B = v1;
        S = 2*(numel(intersect(A,B)))/(numel(A) + numel(B));
        SSS(1,j) = S;

        %solution 2
        A = C(2).vector; %LCBC vtxs
        B = v2;
        S = 2*(numel(intersect(A,B)))/(numel(A) + numel(B));
        SSS(2,j) = S;

        %solution 3
        A = C(3).vector; %LCBC vtxs
        B = v3;
        S = 2*(numel(intersect(A,B)))/(numel(A) + numel(B));
        SSS(3,j) = S;

        %permuted mean across 3 solutions
        pS(j)=mean(SSS(:,j)); 
    end

    %empirical dice
    eD=mm;
    x=eD(c);
    fprintf('\nempirical Dice for %s is = %6.10f',clab,x)
    figure(1)
    histogram(pS)
    xlim([0.09,x+0.05])
    line([x,x],[0,200],'LineWidth',5,'Color','r')

    pdice=mean(pS);
    fprintf('\ntrue expected Dice at random for %s is = %6.10f',clab,pdice)

    %calculate proportion of values that are larger than empirical
    %adding 1 to numerator and denominator
    p=(sum(abs(pS) > abs( eD(c) )) + 1) / (length(pS) + 1);
    fprintf('\np-value for %s = %6.10f',clab,p)
end
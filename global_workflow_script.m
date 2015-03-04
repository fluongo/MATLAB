%initialize distributed computing toolbox
%matlabpool

%Check out a curve fitting toolbox to ensure you have it for calculating F0
 sig = rand(10,1000) 
 temp_ppc = csaps([1:size(sig,2)],sig,0.01);
 C = fnval(temp_ppc,[1:length(sig)]); 


%%

clear QPL_all
clear dir_QPL

%fill in the directories here
dir_QPL = { 'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\WT_HighK_QPL_New\1056_S1_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\WT_HighK_QPL_New\1276.S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\WT_HighK_QPL_New\1956_S1_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\WT_HighK_QPL_New\1956_S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\WT_HighK_QPL_New\1956_S3_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Control\B5_S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Control\B6_S1_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Control\B6_S4_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Control\C3_S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\1630_S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\1630_S3_correct'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\1630_slice1_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\B4_S1_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\B4_S3_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\C1_S3_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\C2_S3_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\C7_S1_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\C7_S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\mouse85_S3_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\mouse85_slice1_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_HighK_QPL\Disc1 DN\mouse85_slice2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_Carb_QPL\Control\B2.9_S1_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_Carb_QPL\Control\B2.9_S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_Carb_QPL\Control\C2.2_S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_Carb_QPL\Control\C2.2_S3_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_Carb_QPL\Disc1 DN\B2.3_S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_Carb_QPL\Disc1 DN\B2.7_S1_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_Carb_QPL\Disc1 DN\B2.7_S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_Carb_QPL\Disc1 DN\C2.5_S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_Carb_QPL\Disc1 DN\C2.9_S2_CORRECT'
    'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\Disc1_Carb_QPL\Disc1 DN\C2.9_S3_CORRECT'
    }

fn_data = sprintf('QPL_data_%s.mat',date);
save(['z:\francisco\Analyzed_Data\' fn_data],'dir_QPL')

for kk = 1:length(dir_QPL)

        tempdir = dir_QPL{kk};
        cd(tempdir);

        % Note dsamp vector first element is the time vector down sample whereas
        % the second element is the time downsample
        m = dir;
        for i = 1:10
            if strfind(m(i).name, '.tif')
                fn = m(i).name
                break
            else
                i = i+1;
            end
        end
        
            
        
        framerate = 10; flims = []; nPCs = 100;
        dsamp = [2 4]; badframes = []; outputdir = [];
        [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA_singleframes(fn, flims, nPCs, dsamp, outputdir ,badframes);

        PCuse  = []; mu = 0.2; nIC = 100; ica_A_guess = []; termtol = [10^-6];  maxrounds = [1000]
        [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, PCuse, mu, nIC, ica_A_guess, termtol, maxrounds);

        smwidth = 0;thresh = 10;arealims = [25 400];plotting = 0;
        [ica_segments, segmentlabel, segcentroid, seg_majoraxis] = CellsortSegmentation(ica_filters, smwidth, thresh, arealims, plotting);

        max_length = 40; [ica_segments, seg_majoraxis, segcentroid] = majoraxisfilter(ica_segments, seg_majoraxis, max_length, segcentroid);
        min_distance = 7;[ica_segments, segcentroid] = distancefilter( ica_segments, segcentroid, min_distance );
        subtractmean = 0;flims = []; cell_sig = CellsortApplyFilter_singleframes(fn, ica_segments, flims, movm, subtractmean);

        % **********Note(**************: Remeber to initialize workers if using the
        % distributed computing toolbox (i.e. parfor), just type 'matlabpool' into
        % the command prompt                          
        [cell_sig_f0 cell_sig_diff cell_sig_f_f0] = ksdensity_normalization_parfor(cell_sig, 2000);
        cell_sig = cell_sig - min(min(cell_sig')) + 1;
        [b,a] = butter(2,0.001);
                    for i = 1:size(cell_sig,1)
                        cell_sig_f0_lowpass(i,:) = filtfilt(b,a,cell_sig(i,:));
                        cell_sig_f_f0_fit(i,:) = (cell_sig(i,:) - cell_sig_f0_lowpass(i,:))./cell_sig_f0_lowpass(i,:);
                    end

        sig = cell_sig_f_f0;
        %sig = cell_sig;

        temp_ppc = csaps([1:size(sig,2)],sig,0.01);
        C = fnval(temp_ppc,[1:length(sig)]); 

        thresh = 2.5;thresh2 = 4;spikewin = 10;mindur = 1;
        spiketimes = [];
        cellid = [];
        spt_beg = [];
        spt_end = [];
        s = [];
        b = [];
        e = [];

                    parfor i = 1:size(sig,1);    
                        [s b e] = detectspikes5(C(i,:)', thresh, thresh2, spikewin, mindur);
                        spiketimes = [spiketimes s];
                        spt_beg = [spt_beg b];
                        spt_end = [spt_end e];
                        cellid = [cellid, i*ones(1,length(s))];
                        i
                    end

    % Excluding cells that have mor than 10 concurrent
    temp = sparse(cellid, spiketimes, 1, size(cell_sig,1),size(cell_sig,2));movt_exclude = find(sum(temp) >10);
    for i = 1:numel(movt_exclude);
        temp = find(spiketimes >= movt_exclude(i) -5  &    spiketimes <= movt_exclude(i) + 5); %excludes 1s around event
        spiketimes(temp) = [];
        cellid(temp) = [];
        spt_beg(temp) = [];
        spt_end(temp) = [];
    end

    clear new_spiketimes new_cellid

    new_spiketimes = [];
    new_cellid = [];
 
    for i = 1:numel(spiketimes)
        new_spiketimes = [spt_beg(i):spiketimes(i) new_spiketimes];
        new_cellid = [cellid(i)*ones(1,spiketimes(i) - spt_beg(i) + 1) new_cellid];
        i;
    end

%Excluding anytime more than 33% of cells are coactive
temp = sparse(new_cellid, new_spiketimes, 1, size(cell_sig,1),size(cell_sig,2));movt_exclude = find(sum(temp) > size(cell_sig,1)/3);
    for i = 1:numel(movt_exclude);
        temp = find(new_spiketimes >= movt_exclude(i) -5  &    new_spiketimes <= movt_exclude(i) + 5); %Exclude 1s around event
        new_spiketimes(temp) = [];
        new_cellid(temp) = [];
    end

% Code for saving data file

temp = pwd; temp = strrep(temp, ' ', '_'); temp = strrep(temp, '\', '_'); 
idx = strfind(temp,'2014');
sub = temp(idx+5:end); 
save('results.mat')

QPL_all{kk} = sparse(new_cellid, new_spiketimes, 1, size(cell_sig,1),size(cell_sig,2));

save(['z:\francisco\Analyzed_Data\' fn_data], 'QPL_all' , '-append')

clearvars -except dir_QPL QPL_all fn_data

end


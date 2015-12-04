%initialize distributed computing toolbox
%matlabpool

% Check out a curve fitting toolbox to ensure you have it for calculating F0
% there are only about 7 of these in the keck center and people in the
% brainard lab are always hogging them....

 sig = rand(10,1000) 
 temp_ppc = csaps([1:size(sig,2)],sig,0.01);
 C = fnval(temp_ppc,[1:length(sig)]); 

%%
%fill in the directories here
dir_Experiments = { 'Z:\Scott\C Drive BU 021315\SCOTT\DATA_TO_ANALYZE\WT_HighK_QPL_New\1056_S1_CORRECT'
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


fn_data = sprintf('Raster_data_%s.mat',date); % Name of mat folder where data will be saved
saveDir = 'z:\francisco\Analyzed_Data\' %Directory in which to save data

save([saveDir fn_data],'dir_Experiments') % saving the directory names

for kk = 1:length(dir_Experiments)

        tempdir = dir_Experiments{kk}; cd(tempdir); % Go to directory of images

        % Code below will simply figure out what the files are named i.e
        % X.tif
        m = dir;
        for i = 1:10
            if strfind(m(i).name, '.tif')
                fn = m(i).name
                break
            else
                i = i+1;
            end
        end
        
        %PCA
        framerate = 10; % Framerate of acquisition
        flims = []; % Leave [] to perform analysis on all frames, otherwise input a range
        nPCs = 100; % Number of PCs to extract
        dsamp = [2 4]; % Downsample factors in [Spatial, Temporal]
        badframes = []; % Frames to exclude, ignore for single-photon
        outputdir = []; % Output directory to save PCA results, by default will save to 'cellsort_preprocessed_data)
        [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA_singleframes(fn, flims, nPCs, dsamp, outputdir ,badframes);
        
        %ICA
        PCuse  = []; % Relevant PCs to use can either be a 2 element array [start, stop] or an array of indices
        mu = 0.2; % Balance between spatial and temporal ICA, 0 is purely temporal, 1 is purely spatial
        nIC = 100; % Number of independent components
        ica_A_guess = []; termtol = [10^-6];  maxrounds = [1000] % Parameters for ICA convergence
        [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, PCuse, mu, nIC, ica_A_guess, termtol, maxrounds);

        % SEGMENTATION
        smwidth = 0; % Standard deviation of guassian blur to apply, 0 applies no blur
        thresh = 10; % Threshold above mean to define as a 'segment'
        arealims = [25 400]; % Limit of size (in pixels) durign segmentation [Low High]
        plotting = 0; % Whether to plot the segmentation outputs or not
        [ica_segments, segmentlabel, segcentroid, seg_majoraxis] = CellsortSegmentation(ica_filters, smwidth, thresh, arealims, plotting);
        
        % FILTERING SEGMENTS
        max_length = 40; [ica_segments, seg_majoraxis, segcentroid] = majoraxisfilter(ica_segments, seg_majoraxis, max_length, segcentroid);
        min_distance = 7;[ica_segments, segcentroid] = distancefilter( ica_segments, segcentroid, min_distance );
        subtractmean = 0;flims = []; cell_sig = CellsortApplyFilter_singleframes(fn, ica_segments, flims, movm, subtractmean);

        % **********Note(**************: Remeber to initialize workers if using the
        % distributed computing toolbox (i.e. parfor), just type 'matlabpool' into
        % the command prompt
        
        
        % KERNEL DENSITY ESTIMATE AND SMOOTHING
        windowSize = 2000; % Length (in samples) of window over which to perform kernel density estimate
        [cell_sig_f0 cell_sig_diff cell_sig_f_f0] = ksdensity_normalization_parfor(cell_sig, windowSize);
        cell_sig = cell_sig - min(min(cell_sig')) + 1; % Set all valus to positive
        [b,a] = butter(2,0.001); % Define butterworth filter coefficients
                    for i = 1:size(cell_sig,1)
                        cell_sig_f0_lowpass(i,:) = filtfilt(b,a,cell_sig(i,:)); % A lowpass filtered version of KS estimate
                        cell_sig_f_f0_fit(i,:) = (cell_sig(i,:) - cell_sig_f0_lowpass(i,:))./cell_sig_f0_lowpass(i,:); % Normalized by low pass KS filtered estime
                    end
        
        % Which signal to use for event detection, I usually use the KSdensity normalized, but
        % other options include the low pass filtered version (cell_sig_f_f0_fit)or the raw
        % signal (cell_sig)
        sig = cell_sig_f_f0; 
        
        % SMOOTHING
        temp_ppc = csaps([1:size(sig,2)],sig,0.01); 
        C = fnval(temp_ppc,[1:length(sig)]); 

        
        % EVENT DETECTION
        thresh = 2.5;% First threshold (to identify putative events)
        thresh2 = 4; % Second threshold (to verify events)
        spikewin = 10; % Window over which to look for large deviations (frames)
        mindur = 1; % Minimum event duration (frames)
        
        % initialize arrays
        spiketimes = []; cellid = []; spt_beg = [];
        spt_end = []; s = []; b = []; e = [];

            % Actual event detection on a cell by cell basis
            parfor i = 1:size(sig,1);    
                [s b e] = detectspikes5(C(i,:)', thresh, thresh2, spikewin, mindur);
                spiketimes = [spiketimes s];
                spt_beg = [spt_beg b];
                spt_end = [spt_end e];
                cellid = [cellid, i*ones(1,length(s))];
                i
            end

    % Excluding cells that have more than 10 perfectly timed spikes (i.e.
    % concurrent single frame) as these are likely to represent movement
    % artifact
    temp = sparse(cellid, spiketimes, 1, size(cell_sig,1),size(cell_sig,2));movt_exclude = find(sum(temp) >10);
    for i = 1:numel(movt_exclude);
        temp = find(spiketimes >= movt_exclude(i) -5  &    spiketimes <= movt_exclude(i) + 5); %excludes 1s around event
        spiketimes(temp) = [];
        cellid(temp) = [];
        spt_beg(temp) = [];
        spt_end(temp) = [];
    end

    % STRETCH FROM EVENT ONSET TO EVENT PEAK
    % Fills in a cell as ative from event onset to event peak
    clear new_spiketimes new_cellid
    new_spiketimes = [];new_cellid = [];
    for i = 1:numel(spiketimes)
        new_spiketimes = [spt_beg(i):spiketimes(i) new_spiketimes];
        new_cellid = [cellid(i)*ones(1,spiketimes(i) - spt_beg(i) + 1) new_cellid];
        i;
    end

    %Excluding anytime more than 33% of cells are coactive because these
    %often represent artifact
    temp = sparse(new_cellid, new_spiketimes, 1, size(cell_sig,1),size(cell_sig,2));movt_exclude = find(sum(temp) > size(cell_sig,1)/3);
        for i = 1:numel(movt_exclude);
            temp = find(new_spiketimes >= movt_exclude(i) -5  &    new_spiketimes <= movt_exclude(i) + 5); %Exclude 1s around event
            new_spiketimes(temp) = [];
            new_cellid(temp) = [];
        end

    % Save all variables for access later
    save('results.mat')
    
    % Cretae and append the raster
    Raster_all{kk} = sparse(new_cellid, new_spiketimes, 1, size(cell_sig,1),size(cell_sig,2));

save(['z:\francisco\Analyzed_Data\' fn_data], 'Raster_all' , '-append')

clearvars -except dir_Experiments Raster_all fn_data

end


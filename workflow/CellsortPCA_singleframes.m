function [mixedsig, mixedfilters, CovEvals, covtrace, movm, ...
    movtm] = CellsortPCA_singleframes(fn, flims, nPCs, dsamp, outputdir, badframes)

% [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn, flims, nPCs, dsamp, outputdir, badframes)
%
% CELLSORT
% Read TIFF movie data and perform singular-value decomposition (SVD)
% dimensional reduction.
% 
%   NOTE: This is the nimmerjahn code except modified to work with single
%   tiff frames as opposed to tiff stacks, fn in this case corresponds to
%   the name of the first frame i.e. img_00000_00.tif
%
% Inputs:
%   fn - movie file name. Must be in TIFF format.
%   flims - 2-element vector specifying the endpoints of the range of
%   frames to be analyzed. If empty, default is to analyze all movie
%   frames.
%   nPCs - number of principal components to be returned
%   dsamp - optional downsampling factor. If scalar, specifies temporal
%   downsampling factor. If two-element vector, entries specify temporal
%   and spatial downsampling, respectively.
%   outputdir - directory in which to store output .mat files
%   badframes - optional list of indices of movie frames to be excluded
%   from analysis
%
% Outputs:
%   mixedsig - N x T matrix of N temporal signal mixtures sampled at T
%   points.
%   mixedfilters - N x X x Y array of N spatial signal mixtures sampled at
%   X x Y spatial points.
%   CovEvals - largest eigenvalues of the covariance matrix
%   covtrace - trace of covariance matrix, corresponding to the sum of all
%   eigenvalues (not just the largest few)
%   movm - average of all movie time frames at each pixel
%   movtm - average of all movie pixels at each time frame, after
%   normalizing each pixel deltaF/F
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

tic
fprintf('-------------- CellsortPCA %s: %s -------------- \n', date, fn)

m = dir;

%-----------------------

% Note, I am calling this function regardless because then if I define a
% flims outside of function, I can just refer to the vals/dir_use vectors
% as indices
[dir_use frame_nt nt_full] = tiff_frames(fn);

    if (nargin<2)||(isempty(flims))
        flims = [1,nt_full];
    end
flims
dsamp_time = dsamp(1);
dsamp_space = dsamp(2);



useframes = setdiff((flims(1):flims(2)), badframes);
% nt = length(1:dsamp_time:length(dir_use));
nt = length(flims(1):dsamp_time:flims(2));

if nargin<3 || isempty(nPCs)
    nPCs = min(150, nt);
end
if nargin<4 || isempty(dsamp)
    dsamp = [1,1];
end
if nargin<5 || isempty(outputdir)
    outputdir = [pwd,'/cellsort_preprocessed_data/'];
end
if nargin<6
    badframes = [];
end
if isempty(dir(outputdir))
    mkdir(pwd, '/cellsort_preprocessed_data/')
end
if outputdir(end)~='/';
    outputdir = [outputdir, '/'];
end

% %Downsampling
% if length(dsamp)==1
%     dsamp_time = dsamp(1);
%     dsamp_space = 1;
% else
%     dsamp_time = dsamp(1);
%     dsamp_space = dsamp(2); % Spatial downsample
% end



[fpath, fname] = fileparts(fn);
    if isempty(badframes)
        fnmat = [outputdir, fname, '_',num2str(flims(1)),',',num2str(flims(2)), '_', date,'.mat'];
    else
        fnmat = [outputdir, fname, '_',num2str(flims(1)),',',num2str(flims(2)),'_selframes_', date,'.mat'];
    end
    
    if ~isempty(dir(fnmat))
        fprintf('CELLSORT: Movie %s already processed;', ...
            fn)
        forceload = input(' Re-load data? [0-no/1-yes] ');
        if isempty(forceload) || forceload==0
            load(fnmat)
            return
        end
    end

fncovmat = [outputdir, fname, '_cov_', num2str(flims(1)), ',', num2str(flims(2)), '_', date,'.mat'];

[pixw,pixh] = size(imread(fn,1));
npix = pixw*pixh;

fprintf('   %d pixels x %d time frames;', npix, nt*dsamp_time)
    if nt<npix
        fprintf(' using temporal covariance matrix.\n')
    else
        fprintf(' using spatial covariance matrix.\n')
    end


%% note amending th code to always use tcov


% Create covariance matrix
    if nt < npix
        [covmat, mov, movm, movtm] = create_tcov(fn, pixw, pixh, useframes, nt, dsamp);
    else
    %Note: This shoul;d be xcov and not tcov    

        [covmat, mov, movm, movtm] = create_tcov(fn, pixw, pixh, useframes, nt, dsamp);
    end
covtrace = trace(covmat) / npix;

movm = reshape(movm, pixw/dsamp_space, pixh/dsamp_space);

% if nt < npix

    % Perform SVD on temporal covariance
    [mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix);

    % Load the other set of principal components
	[mixedfilters] = reload_moviedata(npix, mov, mixedsig, CovEvals);
% 
% else
%     % Perform SVD on spatial components
%     [mixedfilters, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix);
% 
%     % Load the other set of principal components
%     [mixedsig] = reload_moviedata(nt, mov', mixedfilters, CovEvals);
% end

        if (dsamp_space==1)
       %     mixedfilters = reshape(mixedfilters, pixw,pixh,nPCs);
            mixedfilters = reshape(mixedfilters, nPCs,pixw,pixh);
        else
           mixedfilters = reshape(mixedfilters, pixw_dsamp,pixh_dsamp,nPCs);
            
           
            
            % my addition to the code to get the filter out in the right size...
            mixedfilters = imresize(mixedfilters,dsamp_space);
            movm = imresize(movm,dsamp_space);
            %RESHAPING OF THE mixesig array to match output of nPCs x T
            
%                 temp = zeros(size(mixedsig,1),nt*dsamp_time);
%                 for i = 1:size(mixedsig, 1)
%                     temp(i,:) = interp(mixedsig(i,:), dsamp_time^2);
%                 end
%                 mixedsig = temp;
    
        end


firstframe_full = imread(fn,1);
firstframe = firstframe_full;
    if dsamp_space>1
        firstframe = imresize(firstframe, size(mov(:,:,1)),'bilinear');
    end

%------------
% Save the output data
save(fnmat,'mixedfilters','CovEvals','mixedsig', ...
    'movm','movtm','covtrace', 'fn', 'flims', 'nPCs', 'dsamp', 'outputdir', 'badframes')
fprintf(' CellsortPCA: saving data and exiting; ')
toc

    function [covmat, mov, movm, movtm] = create_xcov(fn, pixw, pixh, useframes, nt, dsamp)
        %-----------------------
        % Load movie data to compute the spatial covariance matrix

        npix = pixw*pixh;

        % Downsampling
        if length(dsamp)==1
            dsamp_time = dsamp(1);
            dsamp_space = 1;
        else
            dsamp_time = dsamp(1);
            dsamp_space = dsamp(2); % Spatial downsample
        end

        if (dsamp_space==1)
            mov = zeros(pixw, pixh, nt);
            for jjind=1:length(useframes)
                jj = useframes(jjind);
                mov(:,:,jjind) = imread(fn,jj);
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        else
            [pixw_dsamp,pixh_dsamp] = size(imresize( imread(fn,1), 1/dsamp_space, 'bilinear' ));
            mov = zeros(pixw_dsamp, pixh_dsamp, nt);
            for jjind=1:length(useframes)
                jj = useframes(jjind);
                mov(:,:,jjind) = imresize( imread(fn,jj), 1/dsamp_space, 'bilinear' );
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        end

        fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
        toc
        mov = reshape(mov, npix, nt/dsamp_time);

        % DFoF normalization of each pixel
        movm = mean(mov,2); % Average over time
        movmzero = (movm==0);
        movm(movmzero) = 1;
        mov = mov ./ (movm * ones(1,nt)) - 1; % Compute Delta F/F
        mov(movmzero, :) = 0;

        if dsamp_time>1
            mov = filter(ones(dsamp_time,1)/dsamp_time, 1, mov, [], 2);
            mov = downsample(mov', dsamp_time)';
        end

        movtm = mean(mov,2); % Average over space
        clear movmzeros

        c1 = (mov*mov')/size(mov,2);
        toc
        covmat = c1 - movtm*movtm';
        clear c1
    end

function [covmat, mov, movm, movtm] = create_tcov(fn, pixw, pixh, useframes, nt, dsamp)
        %-----------------------
        % Load movie data to compute the temporal covariance matrix
       % npix = pixw*pixh;

        % Downsampling

            dsamp_time = dsamp(1);
            dsamp_space = dsamp(2); % Spatial downsample

       % if (dsamp_space==1)
            
            [pixw_dsamp,pixh_dsamp] = size(imresize( imread(fn,1), 1/dsamp_space, 'bilinear' ));
            
            npix = pixw_dsamp * pixh_dsamp;
            
            mov = zeros(pixw_dsamp, pixh_dsamp, nt);
            
            
           %  mov = zeros(pixw, pixh, nt);
        
           kk = 1;
                % Adding the downsamplign of time into this step...
            for jjind=flims(1):dsamp_time:flims(2)
%                jj = useframes(jjind);

                    temp_fname = m(dir_use(jjind)).name;
                    
               %  mov(:,:,kk) = imread(temp_fname);
                mov(:,:,kk) = imresize( imread(temp_fname), 1/dsamp_space, 'bilinear' );
                
                kk = kk+1;
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt*dsamp_time)
                    toc
                end

            end
         
         
            
%         else
%             [pixw_dsamp,pixh_dsamp] = size(imresize( imread(fn,1), 1/dsamp_space, 'bilinear' ));
%             mov = zeros(pixw_dsamp, pixh_dsamp, nt);
% 
%             for jjind=1:length(useframes)
%                 jj = useframes(jjind);
%                 mov(:,:,jjind) = imresize( imread(fn,jj), 1/dsamp_space, 'bilinear' );
%                 if mod(jjind,500)==1
%                     fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt*dsamp_time)
%                     toc
%                 end
%  
%                 npix = pixw_dsamp *pixh_dsamp;
%  
%             end
       % end

        fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt*dsamp_time)
        toc
        mov = reshape(mov, npix, nt);

        
        movm = mean(mov,2); % Average over time
        
        
%         % DFoF normalization of each pixel
%         movm = mean(mov,2); % Average over time
%         movmzero = (movm==0); % Avoid dividing by zero
%         movm(movmzero) = 1;
%         mov = mov ./ (movm * ones(1,nt)) - 1;
%         mov(movmzero, :) = 0;

%         if dsamp_time>1
%             mov = filter(ones(dsamp,1)/dsamp, 1, mov, [], 2);
%             mov = downsample(mov', dsamp)';
%         end

%         if dsamp_time>1
%             mov = filter(ones(dsamp_time,1)/dsamp_time, 1, mov, [], 2);
%             mov = downsample(mov', dsamp_time)';
%         end


        c1 = (mov'*mov)/npix;
        movtm = mean(mov,1); % Average over space

        covmat = c1 - movtm'*movtm;
        clear c1
end
%    end
%%
    function [mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix)
        %-----------------------
        % Perform SVD

        covtrace = trace(covmat) / npix;

        opts.disp = 0;
        opts.issym = 'true';
        if nPCs<size(covmat,1)
            [mixedsig, CovEvals] = eigs(covmat, nPCs, 'LM', opts);  % pca_mixedsig are the temporal signals, mixedsig
        else
            [mixedsig, CovEvals] = eig(covmat);
            CovEvals = diag( sort(diag(CovEvals), 'descend'));
            nPCs = size(CovEvals,1);
        end
        CovEvals = diag(CovEvals);
        if nnz(CovEvals<=0)
            nPCs = nPCs - nnz(CovEvals<=0);
            fprintf(['Throwing out ',num2str(nnz(CovEvals<0)),' negative eigenvalues; new # of PCs = ',num2str(nPCs),'. \n']);
            mixedsig = mixedsig(:,CovEvals>0);
            CovEvals = CovEvals(CovEvals>0);
        end

        
      mixedsig = mixedsig' * nt;
      
   %   mixedsig = interp(mixedsig,dsamp_time^2);
      
        CovEvals = CovEvals / npix;

        percentvar = 100*sum(CovEvals)/covtrace;
        fprintf([' First ',num2str(nPCs),' PCs contain ',num2str(percentvar,3),'%% of the variance.\n'])
    end

    function [mixedfilters] = reload_moviedata(npix, mov, mixedsig, CovEvals)
        %-----------------------
        % Re-load movie data
        nPCs = size(mixedsig,1);

        Sinv = inv(diag(CovEvals.^(1/2)));
   
        
%        movtm = interp(mean(mov,1),dsamp_time^2); % Upsample back to original
%         temp = downsample(movtm,dsamp_time^2);
%   movuse = mov - ones(npix,1) * temp;
%         
         movuse = mov - ones(npix,1) * movtm;
         
        mixedfilters = reshape(movuse * mixedsig' * Sinv, npix, nPCs);

        
    end

% This below function can be used to output the number of frames in the
% given movie.....

    function [directory_use vals j] = tiff_frames(fn)
        
        m = dir;
        directory_use = [];
        vals = zeros(numel(m),1);

        for i = 1:numel(m)
            temp = strfind(m(i).name,'_');
            temp_2 = strfind(m(i).name,'.tif');
            if isempty(temp) || isempty(temp_2)
                vals(i) = 0;
            else
              %  vals(i) = str2double(m(i).name(temp(1)+1 : temp_2(1)-1));
                vals(i) = str2double(m(i).name(temp(1)+1 : temp(2)-1));
                directory_use = [directory_use i];
            end
        end
        
        j = max(vals)
        
    end
 end
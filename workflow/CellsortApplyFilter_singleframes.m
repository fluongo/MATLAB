function cell_sig = CellsortApplyFilter_singleframes(fn, ica_segments, flims, movm, subtractmean)
% cell_sig = CellsortApplyFilter(fn, ica_segments, flims, movm, subtractmean)
%
%CellsortApplyFilter
% Read in movie data and output signals corresponding to specified spatial
% filters
%
% Inputs:
%     fn - file name of TIFF movie file
%     ica_segments - nIC x X matrix of ICA spatial filters
%     flims - optional two-element vector of frame limits to be read
%     movm - mean fluorescence image
%     subtractmean - boolean specifying whether or not to subtract the mean
%     fluorescence of each time frame
%
% Outputs:
%     cell_sig - cellular signals
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%
 [dir_use vals nt] = tiff_frames(fn);
if (nargin<3)||isempty(flims)
    flims = [1,nt];
else
    nt = diff(flims)+1;
end
if nargin<5
    subtractmean = 0;
end

[pixw,pixh] = size(imread(fn,1));
if (nargin<4)||isempty(movm)
    movm = ones(pixw,pixh);
else
    movm = double(movm);
end
movm(movm==0) = 1; % Just in case there are black areas in the average image
 k=0;

cell_sig = zeros(size(ica_segments,1), nt);
ica_segments = reshape(ica_segments, [], pixw*pixh);
tic
fprintf('Loading %5g frames for %d ROIs.\n', nt, size(ica_segments,1))
while k<nt
    ntcurr = min(500, nt-k);
    mov = zeros(pixw, pixh, ntcurr);
    for j=1:ntcurr
        temp_fname = m(dir_use(j+k+flims(1)-1)).name;
        movcurr = imread(temp_fname);
        mov(:,:,j) = movcurr;
    end
 %   mov = (mov ./ repmat(movm, [1,1,ntcurr])) - 1; % Normalize by background and subtract mean

    if subtractmean
        % Subtract the mean of each frame
        movtm = mean(mean(mov,1),2);
        mov = mov - repmat(movtm,[pixw,pixh,1]);
    end

    mov = reshape(mov, pixw*pixh, ntcurr);
    cell_sig(:, k+[1:ntcurr]) = ica_segments*mov;

    k=k+ntcurr;
    fprintf('Loaded %3.0f frames; ', k)
    toc
end

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
                vals(i) = str2double(m(i).name(temp(1)+1 : temp(2)-1));
                directory_use = [directory_use i];
            end
        end
        
        j = max(vals)
        
   end
end
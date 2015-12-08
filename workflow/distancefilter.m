function [ ica_segments_f, segcentroid_f, exclude ] = distancefilter( ica_segments, segcentroid, min_distance )
%DISTANCEFILTER Filter to remove cells that are duplicates across multiple
%ICs, by measuring the distance between segment centroids
% Inputs:
%       ica_segments - Output from ica segmentation
%       segcentroid - Corresponding centroids for each segment
%       min_distance - Minimum distance (in pixels) for which to say two segments are
%       the same

%  Outputs:
%       ica_segments_f - Filtered ica_segments, by default will keep the
%       higher ordered of the two segments, this will generally correspond
%       to the higher SNR segment
%       segcentroid_f - Corresponding centroids
%       exclude - Indices of excluded segments


D = pdist(segcentroid);
sq_D = squareform(D);

[I J] = find(sq_D<min_distance & sq_D ~=0);

a = [I,J];
b = [I,J];
b = sort(b,2);
include = [];
exclude = [];

while numel(b) ~= 0 
    temp = b(1,:);
    
    if numel(find(exclude == temp(1))) > 0
        ind = find(b==temp(1));
        ind_2 = find(b==temp(2));
        b([ind ; ind_2]) = temp(1);
        ii = find(diff(b') ==0);
        b(ii,:) = [];
    else    
        ind = find(b==temp(1));
        ind_2 = find(b==temp(2));
        include = [include temp(1)];
        exclude = [exclude temp(2)];
        b([ind ; ind_2]) = temp(1);
        ii = find(diff(b') == 0);
        b(ii,:) = [];
    
    end
end

include = unique(include)
exclude


segcentroid_f = segcentroid;
ica_segments_f = ica_segments;
segcentroid_f(exclude,:) = [];
ica_segments_f(exclude,:,:) = [];


function [ ica_segments_f, seg_majoraxis_f, segcentroid_f ] = majoraxisfilter(ica_segments, seg_majoraxis, max_length, segcentroid)
%MAJORAXISFILTER A filter for dealing with some long skinny segments that
%would sometimes get pulled out, simply throws out segments that have a
%very large major axis

ica_segments_f = ica_segments;
seg_majoraxis_f = seg_majoraxis;
segcentroid_f = segcentroid;

index = find(seg_majoraxis > max_length);

ica_segments_f(index,:,:) = [];
seg_majoraxis_f(index) = [];
segcentroid_f(index, :) = [];


end


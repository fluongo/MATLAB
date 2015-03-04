function [ ica_segments_f, seg_majoraxis_f, segcentroid_f ] = majoraxisfilter(ica_segments, seg_majoraxis, max_length, segcentroid)
%MAJORAXISFILTER Summary of this function goes here
%   Detailed explanation goes here

ica_segments_f = ica_segments;
seg_majoraxis_f = seg_majoraxis;
segcentroid_f = segcentroid;

index = find(seg_majoraxis > max_length);

ica_segments_f(index,:,:) = [];
seg_majoraxis_f(index) = [];
segcentroid_f(index, :) = [];


end


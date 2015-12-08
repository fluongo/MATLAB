function [spiketimes spt_beg spt_end] = detectspikes5(C, thresh, thresh2, spikewin, ...
                                      mindur)

    % C is the calcium trace to be analyzed
    %  [spiketimes spt_beg spt_end] = detectspikes5(C, thresh, thresh2, spikewin, mindur)
    % thresh is the threshold, in factors of the median
    % point-by-point change in C, for detecting candidate events
    % (usually 3)

    % thresh2 is the threshold, in the same units, for the total
    % amplitude of candidate events (usually 15-20)

    % spikewin is the number of points between the detection and
    % the end of an event
    
    % mindur is the minimum duration in points for an event
    
    % RETURNS:
    %               all outputs are as 1xNspikes
    %           spiketimes = the peak of each spike time
    %           spt_beg = Corresponding start time of that event
    %           spt_end = Corresponding end time
    
N = length(C);

Dmask = ones(N-1,1);

spiketimes = [];
spt_beg = [];
spt_end = [];

% identify the points with the largest changes
D = C(2:N) - C(1:N-1);

medianchange = median(abs(D));
mean(abs(D));


for i=1:N-spikewin+1,
    x = C(i:i+spikewin-1);
    y(i) = var(x);
    z(i) = median(abs(x - mean(x)));
end



mean(y);
median(y);

% medianchange = mean(z);
medianchange = median(z);

D = D .* Dmask;
[maxD, locmaxD] = max(D);
locmaxD = locmaxD + 1;        

while maxD && maxD > thresh*medianchange,
    %for j=1:1,
    
    % find the beginning of the spike event
    k = 1;
    while locmaxD-k-1 >= 1 && C(locmaxD-k) > C(locmaxD-k-1),
        k = k+1;
    end
    locbeg = locmaxD-k;
    
    % find the peak of the spike event
    k = 1;
    while locmaxD+k <= N && C(locmaxD+k) >= C(locmaxD+k-1) && locmaxD+k,
        k = k+1;
    end
    k = k-1;
        
    locpeak = locmaxD + k; % location of the peak
    
    if locpeak+spikewin > N,            
        Dmask(locmaxD-1) = 0;
        D = C(2:N) - C(1:N-1);
        D = D .* Dmask;
        [maxD, locmaxD] = max(D);
        locmaxD = locmaxD + 1;
        continue;
    end
    
    % find the end of the spike event
    if sum(C(locpeak+1:locpeak+spikewin) <= C(locbeg)),
        spikeend = min(find(C(locpeak+1:locpeak+spikewin) <= C(locbeg)));
    else
        spikeend = spikewin;
    end
    
    spikeend = spikeend - 1;
    
    locend = locpeak + spikeend;
    
    % check whether the amplitude and duration of the spike exceed
    % the minimum thresholds -- if so, then keep this spike
    if locend - locpeak > mindur %&& C(locpeak) - C(locbeg) > thresh2*medianchange,
        spiketimes = [spiketimes, locpeak];
        spt_beg = [spt_beg, locbeg];
        spt_end = [spt_end, locpeak];
        Dmask(locbeg:(locpeak-1)) = 0;            
        
    else
        % otherwise indicate that no spike is to be added here
        Dmask(locmaxD-1) = 0;            
    end
    D = C(2:N) - C(1:N-1);
    D = D .* Dmask;
    [maxD, locmaxD] = max(D);
    locmaxD = locmaxD + 1;
end

length(spiketimes);

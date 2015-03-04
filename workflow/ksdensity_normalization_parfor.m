function [cell_sig_f0 cell_sig_diff cell_sig_f_f0] = ksdensity_normalization_parfor(cell_sig, win)
% KS density normalization using a sized window as indicated by win
% [cell_sig_f0 cell_sig_diff cell_sig_f_f0] = ksdensity_normalization(cell_sig, win)
% NOTE: This is the version that attempts to use the parallel processing
% toolbox
cell_sig = cell_sig - min(min(cell_sig')) + 1;
cell_sig_f0 = zeros(size(cell_sig));


for i =   1:size(cell_sig,1)
    i
    parfor j = 1:size(cell_sig,2)
        if j < win/2 + 1
            [f xi] = ksdensity(cell_sig(i,1: j+win/2));
            [temp idx] = max(f);
            % cell_sig_f0(i,j) = xi(idx); 
            ksd(j) = xi(idx);
        elseif j > size(cell_sig,2)-(win/2+1)
            [f xi] = ksdensity(cell_sig(i, j - win/2 : end));
            [temp idx] = max(f);
            % cell_sig_f0(i,j) = xi(idx); 
            ksd(j) = xi(idx);
        else    
            [f xi] = ksdensity(cell_sig(i,j - win/2: j+win/2));
            [temp idx] = max(f);
            % cell_sig_f0(i,j) = xi(idx);
            ksd(j) = xi(idx);
        end
    end
    
    cell_sig_f0(i,:) = ksd;
    
end

cell_sig_diff = cell_sig - cell_sig_f0;
cell_sig_f_f0 = (cell_sig_diff)./cell_sig_f0;
% 
% 
% for i = 1:10
%     figure
%     plot(cell_sig(i,:)); hold on; plot(cell_sig_f0(i,:) - median(cell_sig_f0(i,:)),'r'); 
%     hold on; plot(cell_sig_diff(i,:),'g')
%     hold on; plot(200*cell_sig_f_f0(i,:)-30,'c')
%     legend('original signal', 'f0', 'f - f0', ' (f - f0)/f0 ') 
% end
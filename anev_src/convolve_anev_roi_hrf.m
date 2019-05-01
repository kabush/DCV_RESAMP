function [BLDgen] = convolve_anev_roi_hrf(NEVgen,SR,HRF_d)
%
%-Description
% This function takes in a matrix of generated neural events and
% parameters describing the HRF function and convolves the neural
% events and HRF function neural events
%    
%-Parameters
% NEVgen - true neural events generated by the model 
% SR - sample rate of time-series (Hz)
% HRF_d - HRF dispersion (but it's not actually, it's something else)

%Compute number of ROIS
    [N,Bn] = size(NEVgen);
    
    %Generate Kernel at the proper frequency
    kernel = spm_advanced_hrf(1/SR,HRF_d); 
    
    %Calc simulation steps related to simulation time
    Bk = numel(kernel); 
    
    %Allocate Memory to store model brfs
    BLDgen = zeros(N,Bn+Bk-1);
    
    %Convert neural events to indices
    for curr_node = 1:N

        %Superimpose all kernels into one time-series
        for i = 1:numel(NEVgen(curr_node,:))
            BLDgen(curr_node,i:(i+Bk-1)) = NEVgen(curr_node,i)*kernel.' + BLDgen(curr_node,i:(i+Bk-1));
        end
    
    end

    %Trim the excess Bk-1 time-points from result
    BLDgen = BLDgen(:,1:Bn);
    
end
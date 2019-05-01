function [NEVprn,BLDprn,BLDobs] = observe_roi(NEVgen,BLDgen,FG,FO,SNphys,SNscan,ARphys,ARscan,Ttrans,noiseOn,pruneOn,percentOn)
%
%-Description
% This function takes in a matrix of generated neural events and
% convolved neural events and processes these data according to
% simulation parameters to generate simulated observations.
%
%-Inputs
% NEVgen - true neural events 
% BLDgen - true BOLD signal 
% FG - frequency of generation (Hz)
% FO - frequency of observation (Hz)
% SN - signal to noise ratio  
% Ttrans - number of transient initial timepoints
% noiseOn - use of noise model is turned on (true/false)
% pruneOn - data is to be pruned of transient range (true/false)
% normalize - observed BOLD is to be normalized (true/false)
%
%-Outputs
% NEVprn - pruned (removed transients) true neural events
% BLDprn - pruned (removed transients) true BOLD signal with physiology noise
% BLDobs - BLDprn desampled with scan noise
    
%Get dimensions of input
    [N,T] = size(BLDgen);
    
    %Step 1: Remove Transient Timepoints 
    %===================================================================
    if pruneOn
        BLDprn = BLDgen(:,(Ttrans+1):end);
        NEVprn = NEVgen(:,(Ttrans+1):end);
    else
        BLDprn = BLDgen;
        NEVprn = NEVgen;
    end

    %Step 2: Add physiological variation
    %===================================================================
    if noiseOn
        
        for curr_node = 1:N        
            
            %Build physiology-based noise
            phys = zeros(1,numel(BLDprn(curr_node,:)));
            for i = 1:numel(BLDprn(curr_node,:))
                if i>1
                    phys(i) = ARphys*phys(i-1) + randn(1);
                else
                    phys(i) = randn(1);
                end
            end
            
            %Normalize physiology noise
            zphys = zscore(phys);

            %Compute combined signal
            BLDprn(curr_node,:) = BLDprn(curr_node,:)+ (sqrt(mean(BLDprn(curr_node,:).^2))/SNphys)*zphys;

        end
    end
    
    %Step 3: Desample BOLD and Neural Events
    %===================================================================
    desample_rate = FG/FO;
    [N,TS] = size(BLDprn);
    BLDobs = zeros(N,TS/desample_rate);  %%DEBUG
    for curr_node = 1:N
        for i = 1:numel(BLDobs(curr_node,:))
            %Select range of high-definition data to compress
            begin_id = (i-1).*desample_rate+1;
            end_id = i.*desample_rate;            
            BLDobs(curr_node,i) = BLDprn(curr_node,end_id);
        end
    end

    %Step 4: Add scanner noise
    %===================================================================
    if noiseOn

        for curr_node = 1:N
            
            %Build physiology-based noise
            scan = zeros(1,numel(BLDobs(curr_node,:)));
            for i = 1:numel(BLDobs(curr_node,:))
                if i>1
                    scan(i) = ARscan*scan(i-1) + randn(1);
                else
                    scan(i) = randn(1);
                end
            end
            
            %Normalize scanner noise time-series
            zscan = zscore(scan);
            
            %Scale the scanner noise time-series to the proper SNR (amplitude-based)
            meanBLD = mean(BLDobs(curr_node,:).^2);
            
            %Compute combined signal
            BLDobs(curr_node,:) = BLDobs(curr_node,:)+ (sqrt(mean(BLDobs(curr_node,:).^2))/SNscan)*zscan;
            
        end
    
    end

    %Step 5: Convert BOLD to Percent Signal Change
    %===================================================================
    if percentOn
        
       for curr_node = 1:N
           meanCurrent = mean(BLDobs(curr_node,:));
           BLDobs(curr_node,:) = BLDobs(curr_node,:)- meanCurrent;
           BLDobs(curr_node,:) = BLDobs(curr_node,:)/meanCurrent;
       end
   end
    
end    
    

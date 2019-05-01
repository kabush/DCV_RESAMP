function [TP,FP,TN,FN,SENS,SPEC] = compute_confusion_matrix(threshold,NEVrec,NEVgen)
%
%-Description
% This function expands the deconvolved neural signal back 
% into time-series form in which the true neural events
% were generatated
%
%-Inputs
% threshold - value descriminating true and false predictions
% NEVrec - reconstructed true neural events
% NEVgen - true neural events
%
%-Outputs
% TP - True positives
% FP - False positives
% TN - True negatives
% FN - False negatives

    TP = 0;
    FP = 0;
    TN = 0;
    FN = 0;
    
    T_ind = find(NEVgen==1);
    F_ind = find(NEVgen==0);
    
    T_pred_ind = find(NEVrec>threshold);
    F_pred_ind = find(NEVrec<=threshold);
    
    TP = intersect(T_ind,T_pred_ind);
    FP = intersect(F_ind,T_pred_ind);
    FN = intersect(T_ind,F_pred_ind);
    TN = intersect(F_ind,F_pred_ind);

    SENS = numel(TP)/(numel(TP)+numel(FN));
    SPEC = numel(TN)/(numel(TN)+numel(FP));

end

function [ROCx, ROCy] = compute_roc(NEVrec,NEVgen)
%Identify performance
%============================================

NEVrec = NEVrec - min(NEVrec);
NEVrec = NEVrec/max(NEVrec);

ROCx = [];
ROCy = [];

for threshold=0:0.01:1.0
    [TP,FP,TN,FN,SENS,SPEC] = compute_confusion_matrix(threshold,NEVrec,NEVgen);
    ROCx = [ROCx,1-SPEC];
    ROCy = [ROCy,SENS];
end

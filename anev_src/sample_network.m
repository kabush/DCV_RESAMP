function [fc] = sample_network(ANmu,ANvar,ANrot)
%
%-Description
% Here
%
%-Inputs
% ANmu - 
% ANsigma -
%
%-Outputs
% fc - sample of functional connectivity

ac = randn(size(ANmu));
ac = ac.*ANvar;
ac = ANrot*ac';
fc = sigmoid(ac'+ANmu);






    

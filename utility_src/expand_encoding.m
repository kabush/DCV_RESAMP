function NEVrec = expand_encoding(NEVdcv,FG,FO)
%
%-Description
% This function expands the deconvolved neural signal back 
% into a time-series having the same sample frequency as the
% the sequence of true neural events
%
%-Inputs
% NEVdcv - deconvolved observed BOLD signal
% K - size of the HRF function (kernel at observation frequency)
% FG - frequency of true neural event generation (Hz)
% FO - frequency of observation (Hz)
%
%-Outputs
% NEVrec - reconstructed true neural events.  Note the size
%          of NEVrec will not equal NEVgen because the first
%          FG/FO-1 points cannot be reconstructed via interpolation

    dcv_seq = 1:numel(NEVdcv);
    rec_seq = 1:(FO/FG):numel(NEVdcv);
    NEVrec = interp1(dcv_seq,NEVdcv,rec_seq);

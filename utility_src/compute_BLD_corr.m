function [BLD_corr] = compute_BLD_corr(variant,signal)

    BLD_corr = zeros(size(variant,1),1);
    
    for i=1:size(variant,1)

%%        %
%%        variant_adj(i,:) = variant(i,:)-min(variant(i,:));
%%        variant_adj(i,:) = variant_adj(i,:)/max(variant_adj(i,:));
%%
%%        %
%%        signal_adj(i,:) = signal(i,:)-min(signal(i,:));
%%        signal_adj(i,:) = signal_adj(i,:)/max(signal_adj(i,:));

        BLD_corr(i,1) = corr(variant(i,:)',signal(i,:)');
    end
    

end



    

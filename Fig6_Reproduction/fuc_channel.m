function [h,hk] = fuc_channel(M,K,variance_h,variance_hk,seed1)

    if(~exist('seed1','var'))
        seed1 = 'shuffle'; % Generate random numbers by system time
    end
    rng(seed1); % Fixed random number seed

    %% Build Channel
    % Generate WGN model of SOI Channel:
    % All vector signatures are generated as realizations of a complex 
    % random vector with i.i.d. entries having real and imaginary parts 
    % of mean zero and unit variance.
    h = sqrt(variance_h/2) * sqrt(1/2)* (randn(M,1) + 1j * randn(M,1));
    if K == 1
        hk = zeros(M,1);
    elseif K>1
        hk = sqrt(variance_hk/2)* sqrt(1/2)* (randn(M,K-1) + 1j * randn(M,K-1));
    end
end
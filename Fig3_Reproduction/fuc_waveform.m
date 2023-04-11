function [x_n,x_k] = fuc_waveform(N,K,flag,seed2)
    %% For the sake of ease of notation, and without loss of generality, 
    % we assume P=1, so that the magnitude of the power associated with
    % each parameter is modeled within the corresponding signature vector 
    % h. So the variance_xn_dB = variance_xk_dB = 0 dB.

    %% Set Random Seed
    if(~exist('seed2','var'))
        seed2 = 'shuffle'; % Generate random numbers by system time
    end
    rng(seed2); % Fixed random number seed
    
    %% Build Random Wavaform
    if flag == 1 % Generate WGN model 
        % Generate WGN model of SOI Waveform:
        x_n = sqrt(1/2) * (randn(N,1) + 1j * randn(N,1));
        % Generate WGN model of K-1 Interfering Waveform:
        x_k = sqrt(1/2) * (randn(N,K-1) + 1j * randn(N,K-1));
            % per raw in same time, per column belonging to same waveform 
            % but in different time
    elseif flag == 2 % Generate QPSK model
        x = randi([0 1],2*N,K); 
            % The first N entries of each column are real numbers, and the 
            % last N entries are imaginary numbers. The first column is xn,
            % and the following K-1 column is interference signal.
        % Generate QPSK model of SOI Waveform :
        x_n = sqrt(1/2)*((2*x(1:N,1)-1)+1j*(2*x(N+1:2*N,1)-1));
        % Generate QPSK model of K-1 Interfering Waveform :
        x_k = zeros(N,K-1);
        for i = 2:1:K
            x_k(:,i-1) = sqrt(1/2)*(2*x(1:N,i)-1+1j*(2*x(N+1:2*N,i)-1));
        end
    end
end
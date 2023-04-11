%% The reproduction of Figure 2 by formula (15); Frobenius Norm

clear;clc;close all;

M = 4; % number of demension
K = 3; % number of interfering waveform
D = 2; % reduced-rank dimension

%% channels and Interfering channels
SNR = 5; %% in dB
variance_h_dB = SNR; %% in dB
variance_h = 10^(variance_h_dB/10);
variance_hk_dB = SNR; %% in dB
variance_hk = 10^(variance_hk_dB/10);
seed1 = 'shuffle'; % Set random number seed
[h,hk] = fuc_channel(M,K,variance_h,variance_hk,seed1);

%% Ideal R and M_matrix
[R,M_matrix] = fuc_ideal_R(M,K,h,hk);

B1_form = zeros((115-5)/5+1,1);
B2_form = zeros((115-5)/5+1,1);

timemax = 300; %% Set Average Times

for N = 5:5:115 % number of observed samples and estimators
    for time = 1:1:timemax
        fprintf('N= %d time= %d \n',N,time);
        c = M/N;
        
        %% White Gaussian Noise
        variance_noise = 1; %% noise power
        rng('shuffle');
        v_n = sqrt(variance_noise/2) * (randn(M,N) + 1j * randn(M,N)); 
            % a zero-mean temporarily and spatially white Gaussian noise 
            % process with power = 1(variance_noise) per column in same 
            % time
        
        %% Waveform to be estimated and Interfering waveform
        flag = 2; % Control SOI model 1 represents wgn 2 represents QPSK
        seed2 = 'shuffle'; % Set random number seed
        % Seed2 is an optional parameter
        [x_n,x_k] = fuc_waveform(N,K,flag,seed2);
            % x_n :per row in different time
            % x_k :per row in same time, per column belonging to same 
            % waveform but in different time
        
        %% Interfering plus Noise Model
        noise = v_n; %% per column in same time
        for i = 1:1:N
            for j = 1:1:K-1
                noise(:,i) = noise(:,i) + x_k(i,j)*hk(:,j);
                    %% a typical interference-plus-noise vector that is 
                    % linearly described as.
                    % Where x_k(i,j) represents the value of x_k at the jth    
                    % interference signal and the ith time. noise(:, i) 
                    % represents the noise plus interference signals of all
                    % dimensions at the ith moment. hk(:, j) represents the
                    % signature vector of the jth interference signal.
            end
        end
        
        %% Obsevered data
        y_n = h*x_n.' + noise; %% per column in same time
        
        %% Covariance matrix estimation
        R_hat = zeros(M,M);
        for j = 1:1:N
            R_hat = R_hat + y_n(:,j) * y_n(:,j)';
        end
        R_hat = 1/N * R_hat;
        
        %% Cayley-Hamilton theorem
        M_matrix_hat = [h zeros(M,M-1)];
        for i = 1:1:M-1
            M_matrix_hat(:,i+1) = R_hat^(i)*h;
        end
        
        %% loop
        S_D_hat = M_matrix_hat(:,1:D); % Set a D-dimensional subspace with 
            % full rank and it's a Krylov matrix
        v_hat = S_D_hat' * h;
        
        B_hat_down = S_D_hat' * R * S_D_hat;
        B_hat_up = S_D_hat' * R_hat * S_D_hat;      % conventional
        B_proposed = fuc_B_proposed(M,D,c,R_hat,h); % proposed
        
        B1 = B_hat_down-B_proposed; % proposed
        B2 = B_hat_down-B_hat_up; % conventional
        
        B1_form(N/5) = B1_form(N/5)+trace(B1'*B1); % proposed
        B2_form(N/5) = B2_form(N/5)+trace(B2'*B2); % conventional
        
    end
end

S_D = M_matrix(:,1:D);
B = S_D' * R * S_D;
B_form = trace(B'*B);
Normalized_Error_Norm_1 = 1/timemax*B1_form/B_form; % proposed
Normalized_Error_Norm_2 = 1/timemax*B2_form/B_form; % conventional

%% figure
N = 5:5:115;
figure('NumberTitle','off','Name', ...
    'Normalized Error Norm varies with the number of samples in Logarithmic form');
semilogy(N,Normalized_Error_Norm_2,'--r',N,Normalized_Error_Norm_1,'--b');
hold on;grid on;
xlabel('N (number of samples)');ylabel('Normalized Error Norm');
legend('Conventional','Proposed');
%% The reproduction of Figure 2 by formula (22); Frobenius Norm

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
        
        %% Build i.i.d random vectors whose entries have zero mean real and 
        % imaginary parts with variance 1/2 and bounded higher moments.
        variance_U = 1;
        U = sqrt(variance_U/2) * (randn(M,N) + 1j * randn(M,N));
    
        %% Covariance matrix estimation by formula (22)
        R_hat = 1/N * sqrtm(R) * (U * U') * sqrtm(R)';
        
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
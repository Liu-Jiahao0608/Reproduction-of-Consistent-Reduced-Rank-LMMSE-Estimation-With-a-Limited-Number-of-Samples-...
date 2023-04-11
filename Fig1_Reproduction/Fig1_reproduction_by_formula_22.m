%% Normalized MSE of the proposed and conventional estimators of the entries
% The reproduction of Figure 1 by formula (22)

clear;clc;close all;

M = 8; % number of demension
N = 10; % number of observed samples and estimators
K = 6; % number of interfering waveform
D = 6; % reduced-rank dimension
c = M/N;

%% channels and Interfering channels
SNR = 10; %% in dB
variance_h_dB = SNR; %% in dB
variance_h = 10^(variance_h_dB/10);
variance_hk_dB = SNR+5; %% in dB
variance_hk = 10^(variance_hk_dB/10);
seed1 = 4; % Set random number seed
[h,hk] = fuc_channel(M,K,variance_h,variance_hk,seed1);

%% Ideal R and M_matrix
[R,M_matrix] = fuc_ideal_R(M,K,h,hk);
S_D = M_matrix(:,1:D);
B = S_D' * R * S_D;

B1_entry_form = zeros(D,D);
B2_entry_form = zeros(D,D);

timemax = 100;

for time = 1:1:timemax
    fprintf('time= %d \n',time);
    
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
    S_D = M_matrix(:,1:D);
    S_D_hat = M_matrix_hat(:,1:D); % Set a D-dimensional subspace with full 
        % rank and it's a Krylov matrix
    v_hat = S_D_hat' * h;
    
    B_hat_down = S_D_hat' * R * S_D_hat;
    B_hat_up = S_D_hat' * R_hat * S_D_hat;      % conventional
    B_proposed = fuc_B_proposed(M,D,c,R_hat,h); % proposed
    
    B1 = B_hat_down-B_proposed; % proposed
    B2 = B_hat_down-B_hat_up; % conventional
    
    B1_entry_form = B1_entry_form+abs((B1./B)).^2; % proposed
    B2_entry_form = B2_entry_form+abs((B2./B)).^2; % conventional
end 

Normalized_MSE_1 = 1/timemax*B1_entry_form; % proposed
Normalized_MSE_2 = 1/timemax*B2_entry_form; % conventional

%% figure
i = 1:1:D-1;
[I,J] = meshgrid(i);
figure('Name','Normalized MSE in Logarithmic form','NumberTitle','off');
surf(I,J,Normalized_MSE_1(1:D-1,1:D-1),'Marker','.','MarkerSize',15, ...
    'MarkerEdgeColor','k','FaceColor','flat');
hold on;grid on;
surf(I,J,Normalized_MSE_2(1:D-1,1:D-1),'Marker','.','MarkerSize',15, ...
    'MarkerEdgeColor','k','FaceAlpha',0.3,'FaceColor','flat');
set(gca,'ZScale','log');
xlabel('Row i');ylabel('Column j');zlabel('Normalized MSE');
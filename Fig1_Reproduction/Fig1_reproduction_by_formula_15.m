%% Normalized MSE of the proposed and conventional estimators of the entries
% The reproduction of Figure 1 by formula (15)

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
        % x_k :per row in same time, per column belonging to same waveform but
        % in different time
    
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
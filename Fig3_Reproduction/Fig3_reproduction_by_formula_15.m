%% The reproduction of Figure 3 by formula (15);

clear;clc;close all;

M = 16; % number of demension
N = 20; % number of observed samples and estimators
K = 16; % number of interfering waveform
Dmax = 4;

%% channels and Interfering channels
SNR = 10; %% in dB
variance_h_dB = SNR; %% in dB
variance_h = 10^(variance_h_dB/10);
variance_hk_dB = SNR+5; %% in dB
variance_hk = 10^(variance_hk_dB/10);
seed1 = 'shuffle'; % Set random number seed
[h,hk] = fuc_channel(M,K,variance_h,variance_hk,seed1);

%% Initialize Q1 and Q2
Q1_conventional = zeros(Dmax,1);
Q2_conventional = zeros(Dmax,1);
Q1_proposed = zeros(Dmax,1);
Q2_proposed = zeros(Dmax,1);

%% Ideal R and M_matrix
[R,M_matrix] = fuc_ideal_R(M,K,h,hk);

timemax = 100; %% Set Average Times

for D = 1:1:Dmax % reduced-rank dimension
    for time = 1:1:timemax
        fprintf('D= %d, time= %d \n',D,time);
        c = M/N;
        
        %% White Gaussian Noise
        variance_noise = 1; %% noise power
        rng('shuffle');
        v_n = sqrt(variance_noise/2) * (randn(M,N) + 1j * randn(M,N)); 
            % a zero-mean temporarily and spatially white Gaussian noise 
            % process with power = 1(variance_noise) per column in same 
            % time
        
        %% Waveform to be estimated and Interfering waveform
        flag = 1; % Control SOI model 1 represents wgn 2 represents QPSK
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
        
        % B_hat_down
        B_hat_down = S_D_hat'*R*S_D_hat;
        B_hat_down_inv = (B_hat_down)^-1;
        Q1_temp_bound = real(v_hat' * B_hat_down_inv * v_hat);
        Q1_bound = Q1_temp_bound;

        % conventional
        B_hat_up = S_D_hat' * R_hat * S_D_hat;
        B_hat_up_inv = (B_hat_up)^-1;
        omega_D_conventional = B_hat_up_inv * v_hat;
        w_D_conventional = S_D_hat * omega_D_conventional;
        Q1_temp_conventional = real(w_D_conventional' * h);
        Q2_temp_conventional = real(w_D_conventional' * R *  w_D_conventional);
        Q1_conventional(D) = Q1_conventional(D)+abs(Q1_temp_conventional-Q1_bound).^2;
        Q2_conventional(D) = Q2_conventional(D)+abs(Q2_temp_conventional-Q1_bound).^2;

        % proposed
        B_proposed = fuc_B_proposed(M,D,c,R_hat,h);
        B_proposed_inv = (B_proposed)^-1;
        omega_D_proposed = B_proposed_inv * v_hat;
        w_D_proposed = S_D_hat * omega_D_proposed;
        Q1_temp_proposed = real(w_D_proposed' * h);
        Q2_temp_proposed = real(w_D_proposed' * R *  w_D_proposed);
        Q1_proposed(D) = Q1_proposed(D) + abs(Q1_temp_proposed-Q1_bound).^2;
        Q2_proposed(D) = Q2_proposed(D) + abs(Q2_temp_proposed-Q1_bound).^2;
    end
end
MSE_Q1_conventional = Q1_conventional/timemax;
MSE_Q2_conventional = Q2_conventional/timemax;
MSE_Q1_proposed = Q1_proposed/timemax;
MSE_Q2_proposed = Q2_proposed/timemax;


%% figure
D = 1:1:Dmax;
figure('Name','MSE of Q1(w) and Q2(w)','NumberTitle','off');
semilogy(D,MSE_Q2_conventional,'--r^',D,MSE_Q2_proposed,'-b^',D, ...
    MSE_Q1_proposed,'-bo',D,MSE_Q1_conventional,'--ro');
xlim([1,Dmax]);set(gca,'xtick',1:1:Dmax);
hold on;grid on;
xlabel('rank (D)');ylabel('MSE');
legend('Conventional Q2','Proposed Q2','Proposed Q1','Conventional Q1');
%% The reproduction of Figure 5 by formula (22)

clear;clc;close all;

M = 30; % number of demension
N = 20; % number of observed samples and estimators
K = 15; % number of interfering waveform
Dmax = 4;
timemax = 100; %% Set Average Times

%% Initialize Q1 and Q2
Q1 = zeros(Dmax,1);
Q1_bound = zeros(Dmax,1);
Q1_proposed = zeros(Dmax,1);
Q2_proposed = zeros(Dmax,1);
Q1_conventional = zeros(Dmax,1);
Q2_conventional = zeros(Dmax,1);

%% channels and Interfering channels
SNR = 10; %% in dB
variance_h_dB = SNR; %% in dB
variance_h = 10^(variance_h_dB/10);
variance_hk_dB = SNR+5; %% in dB
variance_hk = 10^(variance_hk_dB/10);
seed1 = 'shuffle'; % Set random number seed
[h,hk] = fuc_channel(M,K,variance_h,variance_hk,seed1);

for D = 1:1:Dmax % reduced-rank dimension

    %% Ideal R and M_matrix
    [R,M_matrix] = fuc_ideal_R(M,K,h,hk);
    w_mmse = R^-1*h;
    S_D =  M_matrix(:,1:D);
    v = S_D'*h;
    B = S_D'*R*S_D;
    B_inv = B^-1;
    Q1(D) = real(v'*B_inv*v);

    for time = 1:1:timemax
        fprintf('D= %d, time= %d \n',D,time);
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
        
        % B_hat_down
        B_hat_down = S_D_hat'*R*S_D_hat;
        B_hat_down_inv = (B_hat_down)^-1;
        Q1_temp_bound = real(v_hat' * B_hat_down_inv * v_hat);
        Q1_bound(D) = Q1_bound(D)+Q1_temp_bound;

        % conventional
        B_hat_up = S_D_hat' * R_hat * S_D_hat;
        B_hat_up_inv = (B_hat_up)^-1;
        omega_D_conventional = B_hat_up_inv * v_hat;
        w_D_conventional = S_D_hat * omega_D_conventional;
        Q1_temp_conventional = real(w_D_conventional' * h);
        Q2_temp_conventional = real(w_D_conventional' * R *  w_D_conventional);
        Q1_conventional(D) = Q1_conventional(D)+Q1_temp_conventional;
        Q2_conventional(D) = Q2_conventional(D)+Q2_temp_conventional;

        % proposed
        B_proposed = fuc_B_proposed(M,D,c,R_hat,h);
        B_proposed_inv = (B_proposed)^-1;
        omega_D_proposed = B_proposed_inv * v_hat;
        w_D_proposed = S_D_hat * omega_D_proposed;
        Q1_temp_proposed = real(w_D_proposed' * h);
        Q2_temp_proposed = real(w_D_proposed' * R *  w_D_proposed);
        Q1_proposed(D) = Q1_proposed(D) + Q1_temp_proposed;
        Q2_proposed(D) = Q2_proposed(D) + Q2_temp_proposed;
        
    end
end

%% Compute Q1 and Q2
Q2 = Q1;
Q1_bound = Q1_bound/timemax;
Q2_bound = Q1_bound;
Q1_proposed = Q1_proposed/timemax;
Q2_proposed = Q2_proposed/timemax;
Q1_conventional = Q1_conventional/timemax;
Q2_conventional = Q2_conventional/timemax;

%% Compute SINR
SINR_optimum = (Q2./(abs(Q1).^2)-1).^-1;
SINR_bound = (Q2_bound./(abs(Q1_bound).^2)-1).^-1;
SINR_proposed = (Q2_proposed./(abs(Q1_proposed).^2)-1).^-1;
SINR_conventional = (Q2_conventional./(abs(Q1_conventional).^2)-1).^-1;

SINR_optimum_dB = 10*log(SINR_optimum);
SINR_bound_dB = 10*log(SINR_bound);
SINR_proposed_dB = 10*log(SINR_proposed);
SINR_conventional_dB = 10*log(SINR_conventional);

%% Compute MSE
MSE_optimum = 1-2*real(Q1)+Q2;
MSE_bound = 1-2*real(Q1_bound)+Q2_bound;
MSE_proposed = 1-2*real(Q1_proposed)+Q2_proposed;
MSE_conventional = 1-2*real(Q1_conventional)+Q2_conventional;

%% figure
D = 1:1:Dmax;
figure('Name','Performance tradeoff of proposed and conventional filter estimators versus rank (M =30,N =20)'...
    ,'NumberTitle','off');
subplot(121);plot(D,SINR_optimum_dB,'--.k',D,SINR_bound_dB,'-b*',...
    D,SINR_proposed_dB,'--gx',D,SINR_conventional_dB,'-R+');
xlim([1,Dmax]);set(gca,'xtick',1:1:Dmax);
hold on;grid on;
legend('Optimum','Bound','Proposed','Conventional');
xlabel('rank (D)');ylabel('Output SINR (dB)');
subplot(122);semilogy(D,MSE_optimum,'--.k',D,MSE_bound,'-b*',...
    D,MSE_proposed,'--gx',D,MSE_conventional,'-R+');
xlim([1,Dmax]);set(gca,'xtick',1:1:Dmax);
hold on;grid on;
legend('Optimum','Bound','Proposed','Conventional');
xlabel('rank (D)');ylabel('MSE');


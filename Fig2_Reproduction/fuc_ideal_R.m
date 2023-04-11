function [R,M_matrix] = fuc_ideal_R(M,K,h,hk)
%% Compute the ideal covariance matrix R of y[n]

%% Waveform to be estimated and Interfering waveform Covariance_matrix
R_xn = h*h';
R_xk = 0;
for i = 1:1:K-1
    R_xk = R_xk + hk(:,i)*hk(:,i)';
end

%% White Gaussian Noise Covariance_matrix
variance_vn = 1;
R_vn = variance_vn * eye(M);

%% Interfering plus Noise Model Covariance_matrix
R_nn = R_vn + R_xk;

%% Ideal Covariance_matrix R
R = R_xn + R_nn;

%% Cayley-Hamilton theorem
M_matrix = [h zeros(M,M-1)];
for i = 1:1:M-1
    M_matrix(:,i+1) = R^(i)*h;
end

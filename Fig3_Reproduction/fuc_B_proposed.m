function [B_proposed] = fuc_B_proposed(M,D,c,R_hat,h)

%% solve formula (28)
lamda_hat = eig(R_hat);
lamda_hat = sort(lamda_hat);
%lamda_hat = sort(real(lamda_hat));
sum1 = 0;
syms miu_var
for i = 1:1:M
    sum1 = sum1 + lamda_hat(i)/(lamda_hat(i)-miu_var);
end
miu = double(solve(sum1/M == 1/c,miu_var));
%miu = real(miu);

%% build formula (27)
I_M = eye(M);
B_proposed = zeros(D,D);
for i = 1:1:D
    for j = 1:1:D
        for l = 1:1:M
            for k = 1:1:M
                inv1 = (R_hat-miu(k)*I_M)^-1;
                inv2 = (R_hat-miu(l)*I_M)^-1;
                B_proposed(i,j) = B_proposed(i,j) + (miu(l)^(j-1) * ...
                    miu(k)^(i-1) * h'*inv1*R_hat*inv2*h) /...
                    (c/M*trace(R_hat*(inv1*inv1))*...
                    c/M*trace(R_hat*(inv2*inv2)));
            end
        end
    end
end
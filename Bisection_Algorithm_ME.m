function [L_] = Bisection_Algorithm_ME(n,m,Sigma,D,Rho,kesai)
%BISECTION 此处显示有关此函数的摘要
%   此处显示详细说明
% lambda_1 表示D的最大特征值
% vv_1 是lambda_1的特征向量
[vv_1,lambda_1] = eigs(D,1);
%
d = n + n;
LB = lambda_1 * (1 + sqrt(vv_1' * Sigma * vv_1)/Rho);
UB = lambda_1 * (1 + sqrt(trace(Sigma))/Rho);
% repeat
kk = 0;
%% while 循环很容易出不去
while kk < 100
    gamma = (UB + LB)/2;
    L = gamma^2 * inv(gamma * eye(d) - D) * Sigma * inv(gamma * eye(d) - D);
    h_gamma = Rho^2 - trace(Sigma' * (eye(d) - gamma * inv(gamma*eye(d) - D))^2);
    if h_gamma < 0
        LB = gamma;
    else
        UB = gamma;
    end
    Delta = gamma * (Rho^2 - trace(Sigma)) - trace(L' * D) + gamma^2 * trace((inv(gamma*eye(d)-D))' * Sigma);
    if (h_gamma > 0) && (Delta < kesai)
        break
    end
    kk = kk+1;
end
L_ = L;
end


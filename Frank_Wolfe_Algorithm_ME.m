function [S_ME_star] = Frank_Wolfe_Algorithm_ME(n,m,Sigma_ME,~,Rho,delta,LAMBDA_ME_xx)
%FRANK_WOLFE 此处显示有关此函数的摘要
%   此处显示详细说明
%% 输入参数
% Sigma,u,Rho,delta
% 返回Sigma的特征值和特征向量\

DD = eig(Sigma_ME);
sigma_min = min(DD(DD~=0));% 取到Sigma的非零最小特征值min(a(a~=0))
sigma_max = (Rho + sqrt(trace(Sigma_ME)))^2;
% 令
CC = (2 * sigma_max^4)/sigma_min^3;
%% 迭代初始化
S = Sigma_ME;
k = 0;
while k < 10
    alpha = 2/(k+2);
    % 计算梯度
    D = [LAMBDA_ME_xx, -LAMBDA_ME_xx; -LAMBDA_ME_xx, LAMBDA_ME_xx];
    kesai = alpha * delta * CC;
    %% 用Bisection 算法求解： L
    [L_] = Bisection_Algorithm_ME(n,m,Sigma_ME,D,Rho,kesai);
    % 更新S
    S = S + alpha * (L_ - S);
    k = k + 1;
    

    
    
end
%% 最终输出S_star
S_ME_star = S;
end


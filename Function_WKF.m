function [xee_W,V_W] = Function_WKF(n,m,F,H,Q,xee_W,yy_W,V_W,Rho,delta)
% @ Copyright Zhengyu Feng @ UESTC.
% @ Date 2021.11.15.
% @ Version V_1.0.
%% Wasserstein KF 算法迭代
% prediction
% Form the pseudo-nominal distribution
u = [F;H*F] * xee_W;% 联合分布预测均值
Sigma = [F;H*F] * V_W * [F;H*F]' + [eye(n); H] * Q * [eye(n); H]';% * R联合分布预测协方差
% 用Frank-Wolfe 算法求解
[S_star] = Frank_Wolfe_Algorithm(n,m,Sigma,u,Rho,delta);

%% 输出
S_xx = S_star(1:n,1:n);
S_xy = S_star(1:n,n+1:n+m);
S_yx = S_star(n+1:n+m,1:n);
S_yy = S_star(n+1:n+m,n+1:n+m);
u_x = u(1:n);
u_y = u(n+1:n+m);
% 更新协方差矩阵
V_W =  S_xx -  S_xy * inv(S_yy) * S_yx;
% 更新估计值
xee_W = S_xy * inv(S_yy) * (yy_W - u_y) + u_x;
end
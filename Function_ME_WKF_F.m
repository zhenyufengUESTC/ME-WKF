function [xee_MEW,V_MEW] = Function_ME_WKF_F(n,m,F,H,Q,R,Rho,delta,xee_MEW,yy_MEW,V_MEW,sigma_MEW)
%MC_WKF_F 此处显示有关此函数的摘要
%   此处显示详细说明
% prediction
% Form the pseudo-nominal distribution
LAMBDA_ME_xx = eye(n);
LAMBDA_ME_yy = ones(m);
u_ME = [F;H*F] * xee_MEW;% 联合分布预测均值
Sigma_ME = [F;H*F] * V_MEW * [F;H*F]' + [eye(n); H] * Q * [eye(n); H]';%  * R  联合分布预测协方差
% 用Frank-Wolfe 算法求解
pp = 0;
xee_MEW = u_ME(1:n)+ 0.1;
%         xee_MC_t = xee_MC ;
S_ME_xx = Sigma_ME(1:n,1:n);
S_ME_xy = Sigma_ME(1:n,n+1:n+m);
u_ME_x = u_ME(1:n);
u_ME_y = u_ME(n+1:n+m);
Sigma_ME(n+1:n+m,n+1:n+m) = H*S_ME_xx*H';
S_ = Sigma_ME;
while pp < 3  %norm(xee_MC - xee_MC_t)/norm(xee_MC_t) <= 0.001
    error_ME_xx = (u_ME_x - xee_MEW);% inv(Br)*(yy_ME  - H*xee_ME)  
    for kk = 1:n
        LAMBDA_ME_xx(kk,kk) = exp(-(error_ME_xx(kk)^2)/(2*sigma_MEW^2));%/(eij_x^2)
    end
    G = S_ME_xy*(inv(H*S_ME_xx*H'))';
    S_(1:n,1:n) = S_ME_xx * LAMBDA_ME_xx;
    S_(1:n,n+1:n+n) = Sigma_ME(1:n,n+1:n+m)*G'*LAMBDA_ME_xx;
    S_(n+1:n+n,1:n) = G*Sigma_ME(n+1:n+m,1:n)*LAMBDA_ME_xx;
    S_(n+1:n+n,n+1:n+n) = G*Sigma_ME(n+1:n+m,n+1:n+m)*G'*LAMBDA_ME_xx;
    %
    [S_ME_star] = Frank_Wolfe_Algorithm_ME(n,m,S_,u_ME,Rho,delta,LAMBDA_ME_xx);
    %
    S_ = S_ME_star;
    % 更新协方差矩阵
    S_MEE_xx = S_(1:n,1:n);
    S_ME_xx = S_MEE_xx*inv(LAMBDA_ME_xx);
    Br = chol(R);
    error_MC_yy = inv(Br)*(yy_MEW  - H*xee_MEW);%
    
    for tt = 1:m 
        LAMBDA_ME_yy(tt,tt) = exp(-(error_MC_yy(tt)^2)/(2*sigma_MEW^2));%/
    end
    K_ME = S_ME_xx * H'*inv(H*S_ME_xx*H' + Br*inv(LAMBDA_ME_yy)*Br');%   
    % 更新估计值 
    xee_MEW = K_ME * (yy_MEW  - u_ME_y) + u_ME_x;
% 	[xx1,Pk_1] = MEE_KF_F(F,xee_ME,S_MC_xx,H,yy_ME,R,Q,sigma_MEE);
    
    pp = pp + 1;
end
V_MEW =  (eye(n)-K_ME*H)*S_ME_xx*(eye(n)-K_ME*H)' + K_ME*R*K_ME';
end


function [S_ME_star] = Frank_Wolfe_Algorithm_ME(n,m,Sigma_ME,~,Rho,delta,LAMBDA_ME_xx)
%FRANK_WOLFE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% �������
% Sigma,u,Rho,delta
% ����Sigma������ֵ����������\

DD = eig(Sigma_ME);
sigma_min = min(DD(DD~=0));% ȡ��Sigma�ķ�����С����ֵmin(a(a~=0))
sigma_max = (Rho + sqrt(trace(Sigma_ME)))^2;
% ��
CC = (2 * sigma_max^4)/sigma_min^3;
%% ������ʼ��
S = Sigma_ME;
k = 0;
while k < 10
    alpha = 2/(k+2);
    % �����ݶ�
    D = [LAMBDA_ME_xx, -LAMBDA_ME_xx; -LAMBDA_ME_xx, LAMBDA_ME_xx];
    kesai = alpha * delta * CC;
    %% ��Bisection �㷨��⣺ L
    [L_] = Bisection_Algorithm_ME(n,m,Sigma_ME,D,Rho,kesai);
    % ����S
    S = S + alpha * (L_ - S);
    k = k + 1;
    

    
    
end
%% �������S_star
S_ME_star = S;
end


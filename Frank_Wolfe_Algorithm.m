function [S_star] = Frank_Wolfe_Algorithm(n,m,Sigma,u,Rho,delta)
%FRANK_WOLFE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% �������
% Sigma,u,Rho,delta
% ����Sigma������ֵ����������
DD = eig(Sigma);
sigma_min = min(DD(DD~=0));% ȡ��Sigma�ķ�����С����ֵmin(a(a~=0))
sigma_max = (Rho + sqrt(trace(Sigma)))^2;
% ��
CC = (2 * sigma_max^4)/sigma_min^3;
%% ������ʼ��
S = Sigma;
uu = u;% û����
k = 0;
while k < 10
    alpha = 2/(k+2);
    G = S(1:n,n+1:n+m) * inv(S(n+1:n+m,n+1:n+m));
    % �����ݶ�
    D = [eye(n),-G]' * [eye(n),-G];
    %
    kesai = alpha * delta * CC;
    %% ��Bisection �㷨��⣺ L
    [L_] = Bisection_Algorithm(n,m,Sigma,D,Rho,kesai);
    % ����S
    S = S + alpha * (L_ - S);
    k = k + 1;
end
%% �������S_star
S_star = S;
end


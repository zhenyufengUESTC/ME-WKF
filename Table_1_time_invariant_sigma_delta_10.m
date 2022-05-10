clear all;close all;
%%
% @ Copyright Zhengyu Feng @ UESTC.
% @ Date 2021.11.15.
% @ Version V_1.0.
%% Wasserstein KF �㷨
% ��ʼ������
% KFģ�Ͳ���
n = 2;
m = 1;
N = 300;% ��������
% Э�������

xe = ones(n,N);
xx = ones(n,N);
Rho = 0.1; % Wasserstein radius
delta = 1;% Tolerance

F = [0.9802, 0.0196; 0, 0.9802];%
% A = 0.95* eye(n);
H =  [1,-1];


%% �㷨���� kk ��
for kk = 1:500
    % ״̬��˹����
    v1 = randn(n,N)*0.1;
    % ���ɹ۲������Ǹ�˹����
    % v = randn(n,N)*1;
    q1 = randn(m,N)*0.1;q2 = randn(m,N)*10;
    v = randn(m,N);%%��˹����
    vp = rand(m,N);
    for jj = 1:m
        for tt = 1:N
            if vp(jj,tt) > 0.95
                v(jj,tt) = q2(jj,tt);
            else
                v(jj,tt) = q1(jj,tt);
            end
        end
    end
    Q = (v1 * v1')/N;%% ����֪��״̬����Э���� Q
    R = (v * v')/N;
    
    % WKF �㷨��ʼ������
    V_W = eye(n);% WD KF �㷨
    xe_W(:,1) = xe(:,1);
    % MCC WKF �㷨��ʼ������
    V_MCW = eye(n);% WD KF �㷨
    xe_MC(:,1) = xe(:,1);
    sigma_MCC = 2;
    % MEE WKF �㷨��ʼ������
    V_MEW = eye(n); V_MEW1 = eye(n); V_MEW2 = eye(n); V_MEW4 = eye(n); V_MEW6 = eye(n); V_MEW8 = eye(n); V_MEW10 = eye(n);% WD KF �㷨
    xe_ME(:,1) = xe(:,1);xe_ME1(:,1) = xe(:,1);xe_ME2(:,1) = xe(:,1);xe_ME4(:,1) = xe(:,1);xe_ME6(:,1) = xe(:,1);xe_ME8(:,1) = xe(:,1);xe_ME10(:,1) = xe(:,1);
    sigma_MEE = 1;
    % MCKF�㷨��ʼ������
    xe_MCKF(:,1) = xe(:,1);
    Pk_MC = eye(n);
    % MEKF�㷨��ʼ������
    xe_MEKF(:,1) = xe(:,1);
    Pk_MEE = eye(n);
    % MSE�㷨��ʼ������
    xee2 = xe(:,1);
    xe2 = xe(:,1);
    Pk2 = eye(n);
    %% ״̬�����Ŷ�
    for ii = 2:N
        %         delta_t = -delta + 2 * delta * rand(1);% [-1,1]����ľ��ȷֲ������*delta_t
        %                 delta_t = 0;
        F = [0.9802,  0.0196+0.099*10; 0, 0.9802];%.1196  .98022 *delta_t
        %         H =  [1 + 0.099*1*delta_t,-1];
        % KF ��ʵֵ
        xx(:,ii) = F * xx(:,ii - 1) + v1(:,ii);
        yy(:,ii) = H * xx(:,ii) +  v(:,ii);
        %% MSE KF �㷨
        % ����֪��״̬����Э���� Q
        yy_MSE =  yy(:,ii);
        xee2 = xe2(:,ii-1);
        [xee2,Pk2] = Function_MSE_KF_F(F,H,xee2,yy_MSE,Pk2,Q,R);
        xe2(:,ii) = xee2;
        Err_MSE_KF(kk,ii) = norm(xe2(:,ii) - xx(:,ii));
        
        
        %% Wasserstein KF �㷨����
        % prediction
        yy_W = yy(:,ii);
        xee_W = xe_W(:,ii-1);
        [xee_W,V_W] = Function_WKF(n,m,F,H,Q,xee_W,yy_W,V_W,Rho,delta);
        % ���¹���ֵ
        xe_W(:,ii) = xee_W;
        % �������
        Err_WKF(kk,ii) = norm(xe_W(:,ii) - xx(:,ii));
        
        %% MEE KFѭ������
        yy_MEE = yy(:,ii);
        xe_MEE = xe_MEKF(:,ii-1);
        [xe_MEE,Pk_MEE] = Function_MEE_KF(n,m,F,xe_MEE,Pk_MEE,H,yy_MEE,R,Q,sigma_MEE);
        xe_MEKF(:,ii) = xe_MEE;
        % �������ֵ
        Err_MEE_KF(kk,ii) = norm(xe_MEKF(:,ii) - xx(:,ii));
        
        %% MEE Wasserstein KF �㷨����
        %sigma = 0.5 %
        yy_ME = yy(:,ii);
        xee_ME = xe_ME(:,ii-1);
        [xee_ME,V_MEW] = Function_ME_WKF_F(n,m,F,H,Q,R,Rho,delta,xee_ME,yy_ME,V_MEW,0.5);
        xe_ME(:,ii) = xee_ME;
        % �������
        Err_ME_WKF(kk,ii) = norm(xe_ME(:,ii) - xx(:,ii));
        
        %sigma = 1 %
        yy_ME1 = yy(:,ii);
        xee_ME1 = xe_ME1(:,ii-1);
        [xee_ME1,V_MEW1] = Function_ME_WKF_F(n,m,F,H,Q,R,Rho,delta,xee_ME1,yy_ME1,V_MEW1,1);
        xe_ME1(:,ii) = xee_ME1;
        % �������
        Err_ME_WKF1(kk,ii) = norm(xe_ME1(:,ii) - xx(:,ii));
        
        %sigma = 2 %
        yy_ME2 = yy(:,ii);
        xee_ME2 = xe_ME2(:,ii-1);
        [xee_ME2,V_MEW2] = Function_ME_WKF_F(n,m,F,H,Q,R,Rho,delta,xee_ME2,yy_ME2,V_MEW2,2);
        xe_ME2(:,ii) = xee_ME2;
        % �������
        Err_ME_WKF2(kk,ii) = norm(xe_ME2(:,ii) - xx(:,ii));
        
        %sigma = 4 %
        yy_ME4 = yy(:,ii);
        xee_ME4 = xe_ME4(:,ii-1);
        [xee_ME4,V_MEW4] = Function_ME_WKF_F(n,m,F,H,Q,R,Rho,delta,xee_ME4,yy_ME4,V_MEW4,4);
        xe_ME4(:,ii) = xee_ME4;
        % �������
        Err_ME_WKF4(kk,ii) = norm(xe_ME4(:,ii) - xx(:,ii));
        
        %sigma = 6 %
        yy_ME6 = yy(:,ii);
        xee_ME6 = xe_ME6(:,ii-1);
        [xee_ME6,V_MEW6] = Function_ME_WKF_F(n,m,F,H,Q,R,Rho,delta,xee_ME6,yy_ME6,V_MEW6,6);
        xe_ME6(:,ii) = xee_ME6;
        % �������
        Err_ME_WKF6(kk,ii) = norm(xe_ME6(:,ii) - xx(:,ii));
        
        % sigma = 8 %
        yy_ME8 = yy(:,ii);
        xee_ME8 = xe_ME8(:,ii-1);
        [xee_ME8,V_MEW8] = Function_ME_WKF_F(n,m,F,H,Q,R,Rho,delta,xee_ME8,yy_ME8,V_MEW8,8);
        xe_ME8(:,ii) = xee_ME8;
        % �������
        Err_ME_WKF8(kk,ii) = norm(xe_ME8(:,ii) - xx(:,ii));
        
        % sigma = 10 %
        yy_ME10 = yy(:,ii);
        xee_ME10 = xe_ME10(:,ii-1);
        [xee_ME10,V_MEW10] = Function_ME_WKF_F(n,m,F,H,Q,R,Rho,delta,xee_ME10,yy_ME10,V_MEW10,10);
        xe_ME10(:,ii) = xee_ME10;
        % �������
        Err_ME_WKF10(kk,ii) = norm(xe_ME10(:,ii) - xx(:,ii));
        
        
    end
%     fprintf('%d-th Iteration...\n',kk);
end
save Table_1_time_invariant_sigma_delta_10.mat
%% ��ͼ

% figure; hold on;%
% plot(10*log10(mean(Err_WKF)),'g-+');plot(10*log10(mean(Err_ME_WKF1)),'r-+');
% plot(10*log10(mean(Err_MSE_KF)));plot(10*log10(mean(Err_MEE_KF)));%
% legend('Err-WKF','Err-ME-WKF','Err-MSE-KF','Err-MEE-KF');%

% [mean(mean(Err_WKF)),mean(mean(Err_ME_WKF)),mean(mean(Err_MSE_KF)),mean(mean(Err_MCC_KF)),mean(mean(Err_MEE_KF))]%
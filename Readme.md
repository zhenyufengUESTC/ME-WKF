%% 运行结果
% You can run in your PC as the following example:
% first run the 'Table_run.m' file
run('Table_run.m'); % it will cost long times.

% You can load the ' .mat' files to check the results of the Tables,
% for example:
load(Table_1_time_invariant_sigma_delta_5.mat);  
[mean(mean(Err_MSE_KF)),mean(mean(Err_MEE_KF)),mean(mean(Err_WKF)),mean(mean(Err_ME_WKF)),mean(mean(Err_ME_WKF1)),mean(mean(Err_ME_WKF2)),mean(mean(Err_ME_WKF4)),mean(mean(Err_ME_WKF6)),mean(mean(Err_ME_WKF8)),mean(mean(Err_ME_WKF10))] %#ok<VUNUS>

%% Then plot the MSD error as
figure; hold on;%
plot(10*log10(mean(Err_MSE_KF)));plot(10*log10(mean(Err_MEE_KF)),'r-+');plot(10*log10(mean(Err_WKF)),'g-+');%
plot(10*log10(mean(Err_ME_WKF)));plot(10*log10(mean(Err_ME_WKF1)));plot(10*log10(mean(Err_ME_WKF2)));
plot(10*log10(mean(Err_ME_WKF4)));plot(10*log10(mean(Err_ME_WKF6)));plot(10*log10(mean(Err_ME_WKF8)));
plot(10*log10(mean(Err_ME_WKF10)));

legend('Err-MSE-KF','Err-MEE-KF','Err-WKF','Err_ME_WKF_0.5','Err_ME_WKF1','Err_ME_WKF2','Err_ME_WKF4','Err_ME_WKF6','Err_ME_WKF8','Err_ME_WKF10');%

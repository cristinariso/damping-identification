%% Test file for moving block approach using prescribed signal

% Cristina Riso
% criso@gatech.edu

clearvars 
close all
clc

% add paths from previous folder
addpath('..\00_Utilities');


%% Problem parameters (from Bousman and Winkler, 1981) 

% amplitude
A = 1000.0; 

% undamped frequency (Hz and rad/s)
f_1 = 8.0; omega_1 = f_1*2.0*pi;
f_2 = 6.5; omega_2 = f_2*2.0*pi;

% undamped period (s)
T_1 = 1.0/f_1;
T_2 = 1.0/f_2;

% viscous damping factor
zeta_1 = 0.01; 
zeta_2 = 0.02;

% damped frequency (Hz and rad/s)
omega_d_1 = omega_1*sqrt(1.0-zeta_1^2);
omega_d_2 = omega_2*sqrt(1.0-zeta_2^2);

% damped period (s)
T_d_1 = 2.0*pi/omega_d_1;
T_d_2 = 2.0*pi/omega_d_2;

% time discretization (s)
t_i = 0.0; t_f = 10.0; dt = 0.001; t = t_i:dt:t_f;


%% Compute signal 

% frequency component 1
x_1 = A*exp(-zeta_1*omega_1.*t).*sin(omega_d_1*t);

% frequency component 2
x_2 = A*exp(-zeta_2*omega_2.*t).*sin(omega_d_2*t);

% plot frequency component 1
figure(1); plot(t,x_1,'LineWidth',1);
xlabel('Time (s)'); ylabel('Amplitude (-)');
axis([t_i t_f -A A]); hold all;

% plot frequency component 2
figure(2); plot(t,x_2,'LineWidth',1);
xlabel('Time (s)'); ylabel('Amplitude (-)');
axis([t_i t_f -A A]); hold all;


%% Define moving block parameters

% minimum peak
peak_min = A/100;

% peaks from start and end
peak_from_start = 1; peak_from_end = 0;

% block size (must be a multiple of 2)
block_size = 2^9;

% size ratio bounds
size_ratio_lb = 0.25;
size_ratio_ub = 0.50;

% signal selection
signal_index = 2;


%% Apply moving blocks method to signal manually

% select signal
switch signal_index
    case 1
        x = x_1; T_d = T_d_1; T = T_1; omega_d = omega_d_1; omega = omega_1; zeta = zeta_1;
    case 2
        x = x_2; T_d = T_d_2; T = T_2; omega_d = omega_d_2; omega = omega_2; zeta = zeta_2;
end

% compute reference damping
lambda = -zeta*omega;

% find peaks 
[~,locs] = findpeaks(x,'MinPeakHeight',peak_min); 

% cut signal 
t_test = t(locs(peak_from_start):locs(end-peak_from_end));
x_test = x(locs(peak_from_start):locs(end-peak_from_end));

% plot cut signal
figure(signal_index); plot(t_test,x_test,'--','LineWidth',1);

% sample size
sample_size = length(t_test);

% block to sample size ratio
size_ratio = block_size/sample_size;

% adjust block size if too small
while (size_ratio < size_ratio_lb)

    block_size = block_size*2;
    size_ratio = block_size/sample_size;

end

% adjust block size if too large
while (size_ratio > size_ratio_ub)

    block_size = block_size/2;
    size_ratio = block_size/sample_size;

end

% number of blocks 
n_blocks = sample_size-block_size+1;

% allocate moving block function
lnXbar = zeros(1,n_blocks); 

% loop the blocks
for i = 1:n_blocks

    % isolate signal in current block
    x_i = x_test(i:i+block_size-1);

    % compute FFT
    X_i = fft(x_i);
    
    % compute magnitude of FFT
    Xbar_i = abs(X_i/block_size);
    
    % convert to single-sided magnitude spectrum
    Xbar_i = Xbar_i(1:block_size/2+1);
    Xbar_i(2:end-1) = 2.0*Xbar_i(2:end-1);    

    % store natural logarithm of moving block function
    lnXbar(i) = log(max(Xbar_i));

end

% time interval
tau = t_test(1:n_blocks);

% linearly fit moving block function 
coeff = polyfit(tau,lnXbar,1);

% plot moving block function
figure(11); hold all;
plot(tau,lnXbar,'LineWidth',1);
plot(tau,polyval(coeff,tau),'--','Linewidth',1);
xlabel('Block start time (s)'); ylabel('Log of moving block function (-)');

% obtain damping
lambda_approx = coeff(1);


%% Verification of manually computed values 

% verification of damping
err_lambda = (lambda_approx-lambda)/lambda*100;


%% Apply moving block method to signal using functions

% compute recovery rates using matrix pencil method
[~, ~, lambda_approx_1, ~] = compute_damping_mov_block(t_test,x_test,block_size,size_ratio_lb,size_ratio_ub);

% verification of recovery rate 
err_lambda_1 = (lambda_approx_1-lambda)/lambda*100;


%% Store results in table

% error table
error_table = [lambda 0.0; lambda_approx err_lambda; lambda_approx_1 err_lambda_1];
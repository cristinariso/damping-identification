%% Test file for logarithmic decrement approach using prescribed signal

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


%% Define logarithmic decrement parameters

% minimum peak
peak_min = A/100;

% peaks from start and end
peak_from_start = 1; peak_from_end = 0;

% signal selection
signal_index = 1;


%% Apply logarithmic decrement to signal manually

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
[r_peaks,locs] = findpeaks(x,'MinPeakHeight',peak_min); t_peaks = t(locs);

% select peaks
r_peaks = r_peaks(peak_from_start:end-peak_from_end);
t_peaks = t_peaks(peak_from_start:end-peak_from_end);

% plot peaks
figure(signal_index); 
plot(t_peaks,r_peaks,'ko','MarkerFaceColor','k','MarkerSize',6);

% compute period of damped motion
T_d_approx = t_peaks(end)-t_peaks(end-1);

% compute frequency of damped motion
omega_d_approx = 2.0*pi/T_d_approx;

% extract selected peaks 
r_peak_1 = r_peaks(1);
r_peak_2 = r_peaks(end);

% compute number of cycles between peaks
n_cycles = length(r_peaks)-1;

% compute logarithmic decrement
delta_approx = log(r_peak_1/r_peak_2)/n_cycles;

% compute viscous damping factor
zeta_approx = delta_approx/sqrt(delta_approx^2+4.0*pi^2);

% compute undamped frequency
omega_approx = omega_d_approx/sqrt(1.0-zeta_approx^2);

% compute damping (coefficient of exponential decay)
lambda_approx = -zeta_approx*omega_approx;


%% Verification of manually computed values 

% verification of damping
err_lambda = (lambda_approx-lambda)/lambda*100;


%% Apply logarithmic decrement to signal using functions

% compute peaks for logathmic decrement calculation
[t_peaks, r_peaks] = compute_peaks(t,x,0.0,'one','upper',peak_min,peak_from_start,peak_from_end);

% plot peaks
plot(t_peaks,r_peaks,'rs','MarkerFaceColor','r','MarkerSize',4);

% compute damping using logarithmic decrement
lambda_approx_1 = compute_damping_log_dec(t_peaks,r_peaks,1,0);

% verification of damping 
err_lambda_1 = (lambda_approx_1-lambda)/lambda*100;
%% Test file for logarithmic decrement method

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
f_n_1 = 8.0; omega_n_1 = f_n_1*2.0*pi;
f_n_2 = 6.5; omega_n_2 = f_n_2*2.0*pi;

% undamped period (s)
T_1 = 1.0/f_n_1;
T_2 = 1.0/f_n_2;

% viscous damping factor
zeta_1 = 0.01; 
zeta_2 = 0.02;

% damped frequency (Hz and rad/s)
omega_d_1 = omega_n_1*sqrt(1.0-zeta_1^2);
omega_d_2 = omega_n_2*sqrt(1.0-zeta_2^2);

% damped period (s)
T_d_1 = 2.0*pi/omega_d_1;
T_d_2 = 2.0*pi/omega_d_2;

% time discretization (s)
t_i = 0.0; t_f = 10.0; dt = 0.001; t = t_i:dt:t_f;


%% Compute signal    

% signal selection
signal_index = 2;

% select signal
switch signal_index
    case 1
        x = A*exp(-zeta_1*omega_n_1.*t).*sin(omega_d_1*t); omega_d = omega_d_1; omega_n = omega_n_1; zeta = zeta_1;
    case 2
        x = A*exp(-zeta_2*omega_n_2.*t).*sin(omega_d_2*t); omega_d = omega_d_2; omega_n = omega_n_2; zeta = zeta_2;
end

% compute reference damping (coefficient of exponential decay)
lambda = -zeta*omega_n;


%% Logarithmic decrement method parameters

% minimum peak
peak_min = A/100;

% peaks from start and end
peak_from_start = 1; peak_from_end = 0;


%% Plot signal

% plot selected signal
fig1 = figure(1); hold all; set(fig1,'Position',[0 0 1200 900]); 
plot(t,x,'LineWidth',1);
xlabel('Time (s)','Interpreter','latex'); 
ylabel('Amplitude','Interpreter','latex');
axis([t_i t_f -A A]); ax = gca; ax.FontSize = 32;
f = gcf; exportgraphics(f,'signal.png','Resolution',300);


%% Apply logarithmic decrement to signal manually

% find peaks 
[r_peaks,locs] = findpeaks(x,'MinPeakHeight',peak_min); t_peaks = t(locs);

% select peaks
r_peaks = r_peaks(peak_from_start:end-peak_from_end);
t_peaks = t_peaks(peak_from_start:end-peak_from_end);

% plot peaks
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
zeta_approx_1 = delta_approx/sqrt(delta_approx^2+4.0*pi^2);

% compute undamped frequency
omega_n_approx_1 = omega_d_approx/sqrt(1.0-zeta_approx_1^2);

% compute damping (coefficient of exponential decay)
lambda_approx_1 = -zeta_approx_1*omega_n_approx_1;


%% Verification of manually computed values 

% verification of frequency
err_omega_n = (omega_n_approx_1-omega_n)/omega_n*100.0;

% verification of damping ratio 
err_zeta = (zeta_approx_1-zeta)/zeta*100.0;

% verification of damping (coefficient of exponential decay)
err_lambda = (lambda_approx_1-lambda)/lambda*100;


%% Apply logarithmic decrement to signal using functions

% compute peaks for logathmic decrement calculation
[t_peaks, r_peaks] = compute_peaks(t,x,0.0,'one','upper',peak_min,peak_from_start,peak_from_end);

% plot peaks
plot(t_peaks,r_peaks,'rs','MarkerFaceColor','r','MarkerSize',4);

% compute damping using logarithmic decrement
[omega_n_approx_2, zeta_approx_2, lambda_approx_2] = compute_damping_log_dec(t_peaks,r_peaks,1,0);


%% Verification of values computed using function

% verification of frequency
err_omega_n_2 = (omega_n_approx_2-omega_n)/omega_n*100.0;

% verification of damping ratio 
err_zeta_2 = (zeta_approx_2-zeta)/zeta*100.0;

% verification of damping (coefficient of exponential decay)
err_lambda_2 = (lambda_approx_2-lambda)/lambda*100;
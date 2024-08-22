%% Test file for matrix pencil approach using prescribed signal

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


%% Define matrix pencil method parameters

% minimum peak
peak_min = A/100;

% peaks from start and end
peak_from_start = 1; peak_from_end = 0;

% singular value decomposition tolerance
svd_tol = 0.1;

% window size (-1 to consider entire signal at once) 
window_size = -1;

% stride size
stride_size = -1;

% number of samples between peaks
n_samples = 10;

% KS aggregation across modes
ks_rho = 10000;

% minimum frequency
freq_min = 0.0000001;

% signal selection
signal_index = 1;


%% Apply matrix pencil method to signal manually

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

% cut signal 
t_test = t(locs(peak_from_start):locs(end-peak_from_end));
x_test = x(locs(peak_from_start):locs(end-peak_from_end));

% plot cut signal
figure(signal_index); plot(t_test,x_test,'--','LineWidth',1);

% check flag
if (window_size == -1)

    % set window size to full signal
    window_size = length(t_peaks); 
    
    % set stride size to half the window size
    stride_size = floor(window_size/2);

end

% check number of samples
if (n_samples == -1)

    % use all output data
    spacing = 1;

else

    % get original sampling time step
    dt = t_test(2)-t_test(1);

    % compute resampling time step
    dt_spacing = (t_peaks(end)-t_peaks(end-1))/n_samples;

    % compute spacing between output data points
    spacing = ceil(dt_spacing/dt);

end

% plot peaks
figure(signal_index); 
plot(t_peaks,r_peaks,'ko','MarkerFaceColor','k','MarkerSize',6);

% build the matrices indicating the indices of the various windows
mod1 = repmat([1,window_size],floor(length(t_peaks)/stride_size)-1,1); sz = size(mod1);
mod2 = repmat((0:stride_size:length(t_peaks)-2*stride_size)',1,sz(2));

% indices of various windows
vals = mod1+mod2;

% allocations for recovery rate and amplitudes
r_rates = zeros(sz(1),1); rates = zeros(sz(1),1);

% loop the windows
for i = 1:sz(1)

    % get indices of this window
    indices = vals(i,:);

    % stop when signal is over
    if (indices(end) > length(t_peaks))
        continue
    end

    % find index of window start
    i_1 = find(t == t_peaks(indices(1)));

    % find index of window end
    i_2 = find(t == t_peaks(indices(2)));

    % skip if not found
    if (isempty(i_1) || isempty(i_2))
        continue;
    end

    % extract times for this window
    t_i = t(i_1:spacing:i_2);

    % extract data for this window
    y_i = x(i_1:spacing:i_2);

    % set amplitude to the average of this window
    r_rates(i) = (y_i(1)+y_i(end))/2.0;
    
    % set local damping to the aggregate
    rates(i) = compute(t_i,y_i,ks_rho,svd_tol,freq_min);

end

% extract damping
lambda_approx = rates(end);


%% Verification of manually computed values 

% verification of damping
err_lambda = (lambda_approx-lambda)/lambda*100;


%% Apply matrix pencil method to signal using functions

% compute peaks for logathmic decrement calculation
[t_peaks, r_peaks] = compute_peaks(t,x,0.0,'one','upper',peak_min,peak_from_start,peak_from_end);

% identify index of last peak
index = find(t == t_peaks(end));

% cut signal 
t_test = t(1:index);
x_test = x(1:index);

% plot peaks
plot(t_peaks,r_peaks,'rs','MarkerFaceColor','r','MarkerSize',4);

% group parameters
params.svd_tol = svd_tol;
params.window_size = -1;
params.stride_size = -1;
params.n_samples = n_samples;
params.ks_rho = ks_rho;
params.freq_min = freq_min;

% compute recovery rates using matrix pencil method
[~, rates, ~] = compute_damping_mat_pen(params,t_test,x_test,t_peaks);

% extract damping
lambda_approx_1 = rates(end);

% verification of damping
err_lambda_1 = (lambda_approx_1-lambda)/lambda*100;


%% Helper function for aggregate damping

function c = compute(t, x, ks_rho, svd_tol, freq_min)

% data points below an amplitude threshold must be removed upfront 

% uniform time step
dt = t(2)-t(1);

% number of outputs
N_t = length(t);

% matrix pencil parameter
L = floor(N_t/2)-1;

% allocate Hankel matrix
y = zeros([N_t-L,L+1]);

% populate Hankel matrix
for i = 1:N_t-L
    for j = 1:L+1
        y(i,j) = x(i+j-1); 
    end
end

% conduct SVD
[~,sigma,VT] = svdecon(y);

% singular vectors
VT = VT';

% singular values
sigma = diag(sigma);

% normalize singular values
sigma_norm = sigma/max(sigma);

% singular values to be retained
N_above_tol = length(sigma(sigma_norm > svd_tol));

% number of retained modes (model order)
N_m = min(max(2,N_above_tol),L);

% retained singular vectors
Vhat = VT(1:N_m,:);

% first L rows
V1T = Vhat(:,1:end-1);

% last L rows
V2T = Vhat(:,2:end);

% pseudo inverse
V1inv = pinv(V1T);

% matrix to estimate the eigenvalues and eigenvectors
H = V2T*V1inv;

% eigenvalues of H
[~, lam, ~] = eig(H); lam = diag(lam);

% complex exponents
s = log(lam(1:N_m))/dt;

% real and imaginary parts of complex exponents
alpha = real(s); freqs = imag(s);

% retain real parts for frequencies above threshold
alpha = alpha(freqs > freq_min);

% maximum damping
alpha_max = max(alpha);

% damping aggregate
c = alpha_max+log(sum(exp(ks_rho*(alpha-alpha_max))))./ks_rho;

end


%% Helper function for singular value decomposition

function [U,S,V] = svdecon(X)

% size of input matrix
[m,n] = size(X);

% case of less rows than columns
if  m <= n

    C = X*X';
    [U,D] = eig(C);
    clear C;

    [d,ix] = sort(abs(diag(D)),'descend');
    U = U(:,ix);

    if nargout > 2
        V = X'*U;
        s = sqrt(d);
        V = bsxfun(@(x,c)x./c, V, s');
        S = diag(s);
    end

% other case
else

    C = X'*X;
    [V,D] = eig(C);
    clear C;

    [d,ix] = sort(abs(diag(D)),'descend');
    V = V(:,ix);

    s = sqrt(d);
    S = diag(s);
    U = 0;

end

end
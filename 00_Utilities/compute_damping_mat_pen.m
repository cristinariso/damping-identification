%% Dominant damping using matrix pencil method

% Cristina Riso
% criso@gatech.edu

% based on the function originally written by

% Tarun Golla
% tgolla@gatech.edu


%% Main matrix pencil method function

function [r_rates, rates, params_mat] = compute_damping_mat_pen(params, t, y, t_peaks)

% singular value decomposition tolerance
svd_tol = params.svd_tol;

% number of sliding window peaks
window_size = params.window_size;

% overlap between windows
stride_size = params.stride_size;

% number of samples within dominant oscillation period
n_samples = params.n_samples;

% KS aggregation across modes
ks_rho = params.ks_rho;

% minimum frequency
freq_min = params.freq_min;

% check window size
if (window_size == -1)

    % reset window size to use entire signal
    window_size = length(t_peaks);

    % reset stride size to half the window size
    stride_size = floor(window_size/2);

end

% get original sampling time step
dt = t(2)-t(1);

% check number of samples
if (n_samples == -1)
    
    % compute resampling time step
    dt_spacing = dt;

    % use all output data with no resampling
    spacing = 1;
    
    % compute resampling rate
    resampling_rate = 1;    

else

    % compute resampling time step
    dt_spacing = (t_peaks(end)-t_peaks(end-1))/(n_samples-1);

    % compute spacing between output data points
    spacing = ceil(dt_spacing/dt);

    % compute resampling rate
    resampling_rate = spacing*dt/(t_peaks(end)-t_peaks(end-1));

end

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
    y_i = y(i_1:spacing:i_2);

    % set amplitude to the average of this window
    r_rates(i) = (y_i(1)+y_i(end))/2.0;
    
    % set damping to the aggregate
    rates(i) = compute(t_i,y_i,ks_rho,svd_tol,freq_min);

end

% store parameters used for the calculations
params_mat.svd_tol = svd_tol;
params_mat.window_size = window_size;
params_mat.stride_size = stride_size;
params_mat.n_samples = n_samples;
params_mat.dt = dt;
params_mat.dt_spacing = dt_spacing;
params_mat.dt_resampling = t_i(end)-t_i(end-1);
params_mat.window_duration = t_i(end)-t_i(1);
params_mat.spacing = spacing;
params_mat.resampling_rate = resampling_rate;
params_mat.period = t_peaks(end)-t_peaks(end-1);
params_mat.n_samples_resampling = ceil(params_mat.period/params_mat.dt_resampling+1);
params_mat.ks_rho = ks_rho;
params_mat.freq_min = freq_min;

end


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
%% Dominant damping using moving blocks

% Cristina Riso
% criso@gatech.edu

function [tau, lnXbar, rate, params_mov] = compute_damping_mov_block(t, x, block_size, size_ratio_lb, size_ratio_ub)

% note: the function may need additional testing for signals that have more
% than one frequency component where the other frequencies are close to the
% dominant frequency or their contributions are non-negligible

% sanity check on block size
if (mod(block_size,2) == 1)
    error('Error. \n Block size must be a power of 2.');
end

% sanity check on size ratio lower bound
if (size_ratio_lb < 0.25)
    error('Error. \n Block size must be higher than 0.25.');    
end

% sanity check on size ratio upper bound 
if (size_ratio_ub > 0.5)
    error('Error. \n Block size must be lower than 0.5.');    
end

% sample size
sample_size = length(t);

% sanity check on sample size
if (sample_size < 2)
    error('Error. \n Sample size must be at least 2.');    
end

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

% note: the block size must be a power of 2 for the FFT

% number of blocks 
n_blocks = sample_size-block_size+1;

% allocate natural logarithm of moving-block function
lnXbar = zeros(1,n_blocks); 

% loop the blocks
for i = 1:n_blocks

    % extract current block
    x_i = x(i:i+block_size-1);

    % compute FFT 
    X_i = fft(x_i);
    
    % compute magnitude of FFT
    Xbar_i = abs(X_i/block_size);
    
    % convert to single-sided magnitude spectrum
    Xbar_i = Xbar_i(1:block_size/2+1); Xbar_i(2:end-1) = 2.0*Xbar_i(2:end-1);

    % store natural logaritm of moving-block function
    lnXbar(i) = log(max(Xbar_i));

end

% initial times of blocks
tau = t(1:n_blocks);

% linearly fit natural logarithm of moving-block function
coeff = polyfit(tau,lnXbar,1);

% compute damping (fitting function slope)
rate = coeff(1);

% store parameters used for the calculations
params_mov.block_size = block_size;
params_mov.size_ratio = size_ratio;
params_mov.size_ratio_lb = size_ratio_lb;
params_mov.size_ratio_ub = size_ratio_ub;

end
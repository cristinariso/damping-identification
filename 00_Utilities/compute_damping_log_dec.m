%% Dominant damping using logarithmic decrement

% Cristina Riso
% criso@gatech.edu

function rate = compute_damping_log_dec(t_peaks, r_peaks, peak_from_start, peak_from_end)

% compute period of damped motion
T_d = t_peaks(end-peak_from_end)-t_peaks(end-peak_from_end-1);

% compute frequency of damped motion
omega_d = 2.0*pi/T_d;

% extract selected peaks 
r_peak_1 = r_peaks(peak_from_start);
r_peak_2 = r_peaks(end-peak_from_end);

% compute number of cycles between peaks
n_cycles = length(r_peaks)-peak_from_start-peak_from_end;

% compute logarithmic decrement
delta = log(r_peak_1/r_peak_2)/n_cycles;

% compute viscous damping factor
zeta = delta/sqrt(delta^2+4.0*pi^2);

% compute undamped frequency
omega_n = omega_d/sqrt(1.0-zeta^2);

% compute damping (coefficient of exponential decay)
rate = -zeta*omega_n;

end
%% Signal peaks

% Cristina Riso
% criso@gatech.edu

function [t_peaks, r_peaks] = compute_peaks(t, x, x_e, peak_type, peak_level, peak_min, peak_from_start, peak_from_end)

% note: this function assumes the oscillation is symmetric about x_e and
% should be generalized to handle situations whether that is not the case

% subtract equilibrium
x = x-x_e;

% manipulate signal according to chosen peak options
switch peak_type
    case 'one'
        switch peak_level
            case 'lower'
                x = -x;
        end
    case 'both'
        x = abs(x);
end

% find peaks 
[r_peaks,locs] = findpeaks(x,'MinPeakHeight',peak_min); t_peaks = t(locs);

% select peaks
r_peaks = r_peaks(peak_from_start:end-peak_from_end);
t_peaks = t_peaks(peak_from_start:end-peak_from_end);

% sum equilibrium back
r_peaks = r_peaks+x_e;
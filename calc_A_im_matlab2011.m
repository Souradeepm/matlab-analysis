function A_im = calc_A_im_matlab2011(freq)
% calc_A_im_matlab2011 - Calculate imaginary part of DRT projection matrix
%
% This function builds the imaginary component of the matrix that projects the DRT
% (Distribution of Relaxation Times) into measured impedance values.
% Used in both DRT inversion and forward EIS calculation.
%
% The matrix element A_im(p,q) represents the contribution of relaxation time q
% to the imaginary part of impedance at frequency p.
%
% Mathematical basis:
%   Im(Z) = 0.5*(omega*tau)/(1 + (omega*tau)^2) * log(tau_q+1/tau_q-1)
% where omega = 2*pi*f and tau = 1/f
%
% INPUTS:
%   freq - Frequency vector (Hz) [N x 1]
%
% OUTPUTS:
%   A_im - Imaginary part projection matrix [N x N]

omega = 2 * pi * freq(:);
tau = 1 ./ freq(:);
n = numel(freq);

A_im = zeros(n, n);
for p = 1:n
    for q = 1:n
        if q == 1
            log_term = log(tau(q+1) / tau(q));
        elseif q == n
            log_term = log(tau(q) / tau(q-1));
        else
            log_term = log(tau(q+1) / tau(q-1));
        end
        A_im(p, q) = 0.5 * (omega(p) * tau(q)) / (1 + (omega(p) * tau(q))^2) * log_term;
    end
end
end

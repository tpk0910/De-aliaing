function beamform_image = cpwc_beamform_das(v_interp, t0, t_offset, settings, steer_deg, zero_pad, N_ele, sp)
%CPWC_BEAMFORM_DAS Fast DAS beamforming for a single plane-wave angle.
%
% Timing logic matches PW_compounding.m:
%   tof = delay - t0 + t_offset
%   pos = round(tof * settings.itp_fs) + zero_pad
% and sample v_interp(pos, rx), sum over rx.
%
% Inputs:
%   v_interp : [Nt_itp x Nelem] interpolated channel data (already zero padded)
%   t0       : scalar start time from Field II (seconds)
%   t_offset : scalar group delay compensation (seconds)
%   settings : S_M3_1D settings
%   steer_deg: steering angle in degrees
%   zero_pad : padding length used in z_interp
%
% Output:
%   beamform_image: [Nx x Nz] (same as your original before transettings.element_positionose)
%% Default setting

if nargin < 7 || isempty(N_ele)
        N_ele = settings.N_elements;
end
if nargin < 8 || isempty(sp)
        sp = settings.element_position;
end

Nx   = length(settings.img_x);
Nz   = length(settings.img_z);

x = settings.img_x(:);      % [Nx x 1]
z = settings.img_z(:).';    % [1 x Nz]

theta = steer_deg * pi/180;

% TX distance (meters): x*sin(theta) + z*cos(theta)
tx_dist = x .* sin(theta) + z .* cos(theta);  % [Nx x Nz] (implicit expansion)

beamform_image = zeros(Nx, Nz);

z2 = z.^2;
Nt = size(v_interp, 1);

for rx = 1:N_ele
    xe = sp(rx);

    % RX distance (meters): sqrt((xe - x)^2 + z^2)
    dx2 = (xe - x).^2;     % [Nx x 1]
    rx_dist = sqrt(dx2 + z2); % [Nx x Nz]

    delay = (tx_dist + rx_dist) / settings.c;
    tof   = delay - t0 + t_offset;

    pos = round(tof * settings.itp_fs) + zero_pad;

    % Clamp like your script (<=0 -> 1, >Nt -> 1)
    pos(pos <= 0) = 1;
    pos(pos > Nt) = 1;

    beamform_image = beamform_image + reshape(v_interp(pos(:), rx), Nx, Nz);
end
end

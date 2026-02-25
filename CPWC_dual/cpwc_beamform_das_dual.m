function beamform_image = cpwc_beamform_das_dual(v_interp, t0, t_offset, settings, steer_deg, zero_pad, N_ele_x, sp_x)
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

%% -------------------- Resolve Nx (elements per row) & x positions --------------------
if nargin < 7 || isempty(N_ele_x)
    N_ele_x = settings.N_elements_x;
end

if nargin < 8 || isempty(sp_x)
    sp_x = settings.ele_pos_x;
end
sp_x = sp_x(:).';  % row vector

if numel(sp_x) ~= N_ele_x
    error('cpwc_beamform_das_dual_15D:BadSpX', ...
        'sp_x length (%d) must equal N_ele_x (%d).', numel(sp_x), N_ele_x);
end

%% -------------------- Resolve Ny & y positions --------------------
N_ele_y = 1;
if isfield(settings, 'N_elements_y') && ~isempty(settings.N_elements_y)
   N_ele_y = settings.N_elements_y;
end
if isempty(N_ele_y) || N_ele_y < 1, N_ele_y = 1; end

Nt  = size(v_interp, 1);
Nch = size(v_interp, 2);

% infer Ny if missing but channel count matches
if (N_ele_y == 1) && (Nch ~= N_ele_x) && (mod(Nch, N_ele_x) == 0)
    N_ele_y = Nch / N_ele_x;
end

if Nch ~= N_ele_x * N_ele_y
    error('cpwc_beamform_das_dual_15D:ChannelCountMismatch', ...
        'Expected Nch = N_ele_x*Ny = %d*%d=%d but got %d.', ...
        N_ele_x, N_ele_y, N_ele_x*N_ele_y, Nch);
end

% y slice for x-z image
% - S_DRDP_dual 可能有 img_y 向量（mesh 用），beamforming 需要單一切面 y0
y_img = 0;
if isfield(settings,'img_y0') && ~isempty(settings.img_y0)
    y_img = settings.img_y0;
elseif isfield(settings,'img_y_slice') && ~isempty(settings.img_y_slice)
    y_img = settings.img_y_slice;
elseif isfield(settings,'img_y') && ~isempty(settings.img_y)
    % 如果是標量就用；如果是向量就用 0（避免誤用整個向量）
    if isscalar(settings.img_y)
        y_img = settings.img_y;
    else
        y_img = 0;
    end
end

% element y positions (rows)
sp_y = zeros(1, N_ele_y);
if N_ele_y > 1
    if isfield(settings, 'ele_pos_y') && ~isempty(settings.ele_pos_y)
        sp_y = settings.ele_pos_y(:).';
        if numel(sp_y) ~= N_ele_y
            error('cpwc_beamform_das_dual_15D:BadSpY', ...
                'settings.ele_pos_y length (%d) must equal Ny (%d).', numel(sp_y), N_ele_y);
        end
    elseif isfield(settings, 'element_position_y') && ~isempty(settings.element_position_y)
        sp_y = settings.element_position_y(:).';
        if numel(sp_y) ~= N_ele_y
            error('cpwc_beamform_das_dual_15D:BadSpY', ...
                'settings.element_position_y length (%d) must equal Ny (%d).', numel(sp_y), N_ele_y);
        end
    else
        % auto-generate centered y using pitch_y (preferred)
        if isfield(settings,'pitch_y') && ~isempty(settings.pitch_y)
            pitch_y = settings.pitch_y;
        else
            % fallback: height_y + kerf_y (if available), else use settings.height + settings.kerf
            kerf_y = [];
            height_y = [];
            if isfield(settings,'kerf_y') && ~isempty(settings.kerf_y), kerf_y = settings.kerf_y; end
            if isfield(settings,'height_y') && ~isempty(settings.height_y), height_y = settings.height_y; end
            if isempty(kerf_y)
                if isfield(settings,'kerf') && ~isempty(settings.kerf), kerf_y = settings.kerf; else, kerf_y = 0; end
            end
            if isempty(height_y)
                if isfield(settings,'height') && ~isempty(settings.height), height_y = settings.height; else, height_y = 0; end
            end
            pitch_y = height_y + kerf_y;
        end
        sp_y = ((0:N_ele_y-1) - (N_ele_y-1)/2) * pitch_y;
    end
end

%% -------------------- Imaging grid --------------------
Nx_img = numel(settings.img_x);
Nz_img = numel(settings.img_z);

x = settings.img_x(:);      % [Nx_img x 1]
z = settings.img_z(:).';    % [1 x Nz_img]

theta = steer_deg * pi/180;

% TX distance (plane wave), steering only in x-z
tx_dist = x .* sin(theta) + z .* cos(theta);  % [Nx_img x Nz_img]

z2 = z.^2;
beamform_image = zeros(Nx_img, Nz_img);

%% -------------------- Main sum (channel loop) --------------------
for ch = 1:Nch
    % channel ordering: x runs fastest, then y-row blocks
    iy = ceil(ch / N_ele_x);
    ix = ch - (iy-1) * N_ele_x;

    xe = sp_x(ix);
    ye = sp_y(iy);

    dx2 = (xe - x).^2;            % [Nx_img x 1]
    dy2 = (ye - y_img).^2;        % scalar
    rx_dist = sqrt(dx2 + dy2 + z2); % [Nx_img x Nz_img]

    delay = (tx_dist + rx_dist) / settings.c;
    tof   = delay - t0 + t_offset;

    pos = round(tof * settings.itp_fs) + zero_pad;

    % clamp (match your original behavior: out-of-range -> 1)
    pos(pos <= 0) = 1;
    pos(pos > Nt) = 1;

    beamform_image = beamform_image + reshape(v_interp(pos(:), ch), Nx_img, Nz_img);
end

end
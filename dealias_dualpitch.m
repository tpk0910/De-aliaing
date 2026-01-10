function out = dealias_dualpitch(env1, env2, x_im, z_im, opts)
% dealias_dualpitch.m
% Dual-pitch de-aliasing by per-depth normalization + consistency fusion.
%
% Inputs:
%   env1, env2 : envelope images, size [Nz x Nx] (same size)
%   x_im       : 1xNx lateral axis (meters)
%   z_im       : 1xNz axial axis (meters)
%   opts       : struct (optional)
%
% Output:
%   out.score      : [Nz x Nx] de-aliased score map (0..1)
%   out.norm1/norm2: [Nz x Nx] normalized maps
%   out.xhat_z     : [Nz x 1] estimated x per depth (meters)
%   out.ihat_z     : [Nz x 1] peak value per depth
%   out.global_peak: struct with global max location

    if nargin < 5, opts = struct(); end
    if ~isfield(opts, 'use_exclude_center'), opts.use_exclude_center = true; end
    if ~isfield(opts, 'exclude_center_mm'),  opts.exclude_center_mm  = 1.0;  end % exclude +/-1mm around x=0
    if ~isfield(opts, 'pow'),               opts.pow = 1; end % score = (n1*n2)^pow
    if ~isfield(opts, 'eps'),               opts.eps = 1e-12; end

    assert(all(size(env1)==size(env2)), 'env1/env2 must be same size');

    [Nz, Nx] = size(env1);
    x_im = x_im(:).';  % 1xNx
    z_im = z_im(:).';  % 1xNz

    % ----- per-depth robust normalization -----
    norm1 = zeros(Nz, Nx);
    norm2 = zeros(Nz, Nx);

    for iz = 1:Nz
        p1 = env1(iz, :);
        p2 = env2(iz, :);

        % robust baseline removal
        b1 = median(p1);
        b2 = median(p2);

        p1 = p1 - b1;
        p2 = p2 - b2;

        % clamp negatives
        p1(p1 < 0) = 0;
        p2(p2 < 0) = 0;

        % scale to [0,1]
        s1 = max(p1) + opts.eps;
        s2 = max(p2) + opts.eps;

        n1 = p1 / s1;
        n2 = p2 / s2;

        % optional: exclude center region (to avoid mainlobe dominating)
        if opts.use_exclude_center
            xc = opts.exclude_center_mm * 1e-3; % meters
            mask_center = abs(x_im) <= xc;
            n1(mask_center) = 0;
            n2(mask_center) = 0;
        end

        norm1(iz,:) = n1;
        norm2(iz,:) = n2;
    end

    % ----- consistency fusion -----
    score = (norm1 .* norm2);
    if opts.pow ~= 1
        score = score.^opts.pow;
    end

    % ----- per-depth peak (xhat(z)) -----
    [ihat_z, ixhat] = max(score, [], 2);
    xhat_z = x_im(ixhat).';

    % ----- global peak -----
    [global_val, lin] = max(score(:));
    [izg, ixg] = ind2sub(size(score), lin);

    out = struct();
    out.norm1 = norm1;
    out.norm2 = norm2;
    out.score = score;
    out.xhat_z = xhat_z;
    out.ihat_z = ihat_z;
    out.global_peak = struct( ...
        'x', x_im(ixg), ...
        'z', z_im(izg), ...
        'value', global_val, ...
        'ix', ixg, ...
        'iz', izg ...
    );
end

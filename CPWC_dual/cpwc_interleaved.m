% CPWC_DENSE(2*N_sparse+1)x2_TO_WIDE(N_sparse)x2_INTERLEAVED_NO_DEALIAS_15D
% ------------------------------------------------------------------------
% Goal:
%   Use a dense (N_dense x 2) 1.5D array to synthesize the RF of a
%   target wide-element sparse (N_sparse x 2) array with half-pitch interleaving along x,
%   then reconstruct onto Nyquist x-grid and beamform (no de-alias).
%
% Key idea:
%   Dense RF samples approximate v(t,x) on a finer x-grid (pitch_dense = pitch_sparse/2).
%   To emulate WIDE elements and half-pitch interleaving:
%     Row1: discard the 1st dense element -> (N_dense-1)=2*N_sparse samples, then sum pairs -> N_sparse wide elements
%     Row2: discard the last dense element -> (N_dense-1)=2*N_sparse samples, then sum pairs -> N_sparse wide elements
%   This creates two N_sparse-element rows whose effective centers are half-pitch shifted.
%
%   Merge -> sp_merge (2*N_sparse positions), pitch_merge = pitch_sparse/2
%   Map back to Nyquist grid (sp_nyq) by DFT/IDFT (no de-alias),
%   then CPWC DAS beamforming (dual/1.5D).
% ------------------------------------------------------------------------
clc; clear; close all;

%% ---------- Settings ----------
settings        = S_DRDP_dual;
compound_deg    = -3:1:3;

point_position  = [0 0 30e-3];
point_amplitude = 1;

zero_pad        = 500;
DR              = settings.dynamic_range;

% --- target sparse array size in x (VARIABLE) ---
N_sparse = 32;                 % <<< CHANGE THIS (e.g., 64, 22, 16, 48...)
assert(mod(N_sparse,1)==0 && N_sparse>=2, 'N_sparse must be an integer >= 2');

% Dense must be odd = 2*N_sparse+1 for the (drop-1 then pair) construction
N_dense  = 2*N_sparse + 1;     % dense elements per row in x

% --- choose sparse pitch in x (common: preserve aperture) ---
pitch_sparse = 2*settings.aperture_x / N_sparse;
% If you want specifically e.g. 0.67lambda, replace by:
% pitch_sparse = 0.67 * settings.lambda;

kerf_x = settings.kerf_x;

% Dense pitch is half of sparse pitch -> no-alias ground truth
pitch_dense = pitch_sparse/2;

% NOTE:
% width is only used by Field II geometry in cpwc_one_angle_dual.
% Dense elements should be narrower. Wide-element behavior is emulated by pair-sum in RF domain.
width_dense = pitch_dense - kerf_x;

%% ---------- Field II init ----------
field_init(0);
set_sampling(settings.fs);
set_field('c',  settings.c);
set_field('fs', settings.fs);

%% ---------- Allocate compound image ----------
Nx_img = numel(settings.img_x);
Nz_img = numel(settings.img_z);
compound_merge = zeros(Nx_img, Nz_img);

%% ---------- Build dense x positions (N_dense) ----------
x_dense = make_symmetric_positions(N_dense, pitch_dense, 0);

%% ---------- Build WIDE (N_sparse)x2 interleaved positions ----------
% For pair-sum, the effective wide-element centers are means of two dense centers
[x_row1, x_row2] = make_wide_interleaved_positions_from_dense(x_dense);

% merged interleaved lattice (2*N_sparse points)
sp_merge = sort([x_row1(:); x_row2(:)]).';     % 1 x (2*N_sparse)
pitch_merge = pitch_sparse/2;                  % half-pitch interleaved lattice

% precompute mapping indices ONCE (faster + stable for any N_sparse)
[tf1, idx1] = ismember(x_row1(:), sp_merge(:));
[tf2, idx2] = ismember(x_row2(:), sp_merge(:));
assert(all(tf1) && all(tf2), 'Mapping to sp_merge failed. Check position generation.');
idx1 = idx1(:).';   % row vector
idx2 = idx2(:).';

%% ---------- Nyquist grid (target for beamformer) ----------
sp_nyq  = settings.ele_pos_x;
ks_nyq  = 1/settings.pitch_x;
k_range = ks_nyq/2;

%% ---------- CPWC loop ----------
for k = 1:numel(compound_deg)
    steer_deg = compound_deg(k);

    % 1) Dense acquisition (N_dense x 2) by Field II
    [v_dense, t_dense] = cpwc_one_angle_dual( ...
        settings, steer_deg, point_position, point_amplitude, ...
        N_dense, width_dense, x_dense);

    % 2) Interpolation in time (z_interp)
    v_dense_itp = z_interp(v_dense, settings.fs, settings.itp_ratio, zero_pad);
    Nt = size(v_dense_itp, 1);

    % output RF on Nyquist grid (shape matches your beamformer):
    % [Nt, N_elements_x * N_elements_y]
    v_merge_nyq = zeros(Nt, settings.N_elements_x * settings.N_elements_y);

    % 3) Per time sample, per y-row:
    %    dense (N_dense) -> wide row1 (N_sparse) + wide row2 (N_sparse) -> merge (2*N_sparse) -> DFT/IDFT to Nyquist
    for row = 1:Nt
        for iy = 1:settings.N_elements_y

            % Dense channel block for this y-row (assumes v_dense is [Nt, N_dense*Ny])
            cD = (iy-1)*N_dense + (1:N_dense);

            % Extract dense RF vs x for this y-row
            sig_dense = v_dense_itp(row, cD);  % 1 x N_dense
            assert(mod(numel(sig_dense),2)==1, 'Expected odd dense length (2*N_sparse+1).');

            % --- Pair-sum wide-element emulation ---
            % Row1: drop first -> (N_dense-1)=2*N_sparse, then sum pairs -> N_sparse
            sig2N_r1 = sig_dense(2:end);
            sig_r1   = sig2N_r1(1:2:end) + sig2N_r1(2:2:end);

            % Row2: drop last -> (N_dense-1)=2*N_sparse, then sum pairs -> N_sparse
            sig2N_r2 = sig_dense(1:end-1);
            sig_r2   = sig2N_r2(1:2:end) + sig2N_r2(2:2:end);

            % --- Merge onto interleaved lattice (2*N_sparse points) ---
            sig_merge = zeros(1, numel(sp_merge));
            sig_merge(idx1) = sig_r1;
            sig_merge(idx2) = sig_r2;

            % --- Map merged lattice -> Nyquist grid (no de-alias) ---
            [Ym, fk] = DFT(sig_merge, pitch_merge, sp_merge, k_range);
            cF = (iy-1)*settings.N_elements_x + (1:settings.N_elements_x);
            v_merge_nyq(row, cF) = ks_nyq * IDFT(Ym, fk, sp_nyq, ks_nyq);
        end
    end

    clear v_dense v_dense_itp;

    % 4) Group delay compensation
    t_offset = round((length(settings.excitation) + ...
                      length(settings.impulse_response) + ...
                      length(settings.impulse_response) - 2)/2) / settings.fs;

    % 5) Beamforming
    bf_merge = cpwc_beamform_das_dual(v_merge_nyq, t_dense, t_offset, ...
                                     settings, steer_deg, zero_pad, ...
                                     settings.N_elements_x, sp_nyq);

    compound_merge = compound_merge + bf_merge;
    clear v_merge_nyq bf_merge;
end

%% ---------- Image process & plot ----------
log_image = img_process(compound_merge', DR);

figure;
image(settings.img_x*1e3, settings.img_z*1e3, log_image);
axis image; colormap(gray(DR)); colorbar;
title(sprintf('%dx2 -> %dx2 interleaved', ...
    N_dense, N_sparse));
xlabel('x [mm]'); ylabel('z [mm]');
saveas(gcf, sprintf('%d+%d_CPWC.fig', N_sparse, N_sparse));

%   Beam pattern (lateral)
figure;
beam_pattern = max(log_image, [], 1);
beam_pattern = beam_pattern - max(beam_pattern);
plot(settings.img_x*1e3, beam_pattern);
title("Beam pattern");
axis tight;
saveas(gcf, sprintf('%d+%d_Beam_pattern.fig', N_sparse, N_sparse));

%   Beam pattern (axial)
figure;
beam_pattern_z = max(log_image, [], 2);
beam_pattern_z = beam_pattern_z - max(beam_pattern_z);
plot(settings.img_z*1e3, beam_pattern_z);
title("Beam pattern (Axial)");
axis tight;
saveas(gcf, sprintf('%d+%d_Beam_pattern_(Axial).fig', N_sparse, N_sparse));

%% ---------- helper functions ----------
function x = make_symmetric_positions(N, pitch, x_center)
% positions: (-(N/2)+0.5)*pitch : pitch : ((N/2)-0.5)*pitch
if nargin < 3, x_center = 0; end
idx = (-(N/2)+0.5) : 1 : ((N/2)-0.5);
x = x_center + idx * pitch;
end

function [x_row1, x_row2] = make_wide_interleaved_positions_from_dense(x_dense)
% x_dense: 1 x N_dense, where N_dense = 2*N_sparse + 1 (odd)
% Row1: drop first -> 2*N_sparse, then pairwise average centers -> N_sparse
% Row2: drop last  -> 2*N_sparse, then pairwise average centers -> N_sparse

Nd = numel(x_dense);
assert(mod(Nd,2)==1, 'x_dense length must be odd (Nd = 2*N_sparse+1).');

x2N_r1 = x_dense(2:end);      % length Nd-1 (even)
x2N_r2 = x_dense(1:end-1);    % length Nd-1 (even)

assert(mod(numel(x2N_r1),2)==0, 'After dropping one element, length must be even.');

x_row1 = 0.5 * (x2N_r1(1:2:end) + x2N_r1(2:2:end));
x_row2 = 0.5 * (x2N_r2(1:2:end) + x2N_r2(2:2:end));
end
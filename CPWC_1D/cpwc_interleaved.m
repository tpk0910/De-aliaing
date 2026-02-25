% CPWC_2ROW_INTERLEAVED_NO_DEALIAS_SHOW3 (FIXED)
% ------------------------------------------------------------
% Two rows with same N, interleaved by 0.5*pitch (Delta = pitch_row/2).
% NO de-aliasing. Just interleave two rows to form denser lattice.
%
% Show three images:
%   (1) Row 1 only (raw)
%   (2) Row 2 only (raw)
%   (3) Interleaved merge (no de-alias)
%
% FIX:
%   settings.element_position is read-only -> do NOT overwrite settings.
%   Instead, reconstruct merged lattice back onto Nyquist grid by DFT/IDFT.
% ------------------------------------------------------------
clc; clear; close all;

%% ---------- Settings ----------
settings        = S_DRDP_1D;
compound_deg    = -3:1:3;
point_position  = [0 0 30e-3];
point_amplitude = 1;
zero_pad        = 500;
DR              = settings.dynamic_range;

ratio_row = 1/2;  % 1/2 -> ~32 per row, 1/3 -> ~22 per row (rounding depends)

%% ---------- Field II init ----------
field_init(0);
set_sampling(settings.fs);
set_field('c',  settings.c);
set_field('fs', settings.fs);

%% ---------- Allocate compound images ----------
Nx = numel(settings.img_x);
Nz = numel(settings.img_z);

compound_row1  = zeros(Nx, Nz);
compound_row2  = zeros(Nx, Nz);
compound_merge = zeros(Nx, Nz);

%% ---------- CPWC loop ----------
for k = 1:numel(compound_deg)

    steer_deg = compound_deg(k);

    %% ----- Build SAME-N transducer geometry for two interleaved rows -----
    [N_row, width_row, pitch_row, sp_center] = transducer_sym_width( ...
        settings.N_elements, ratio_row, settings.kerf, settings.aperture_size);

    % Interleaved positions (Delta = pitch_row/2), keep overall centered
    sp1 = sp_center - pitch_row/4;
    sp2 = sp_center + pitch_row/4;

    c = mean([sp1(:); sp2(:)]);
    sp1 = sp1 - c;
    sp2 = sp2 - c;

    %% ----- Acquire RF for each row -----
    [v1, t1] = cpwc_one_angle(settings, steer_deg, point_position, point_amplitude, ...
                             N_row, width_row, sp1);
    [v2, t2] = cpwc_one_angle(settings, steer_deg, point_position, point_amplitude, ...
                             N_row, width_row, sp2);

    % Safety: match lengths
    if size(v1,1) < size(v2,1)
        v1 = [v1; zeros(size(v2,1)-size(v1,1), size(v1,2))];
    elseif size(v2,1) < size(v1,1)
        v2 = [v2; zeros(size(v1,1)-size(v2,1), size(v2,2))];
    end

    %% ----- Time interpolation -----
    v1i = z_interp(v1, settings.fs, settings.itp_ratio, zero_pad);
    v2i = z_interp(v2, settings.fs, settings.itp_ratio, zero_pad);

    Nt = size(v1i,1);

    %% ----- Nyquist grid (target for beamformer) -----
    sp_nyq  = settings.element_position;
    ks_nyq  = 1/settings.pitch;
    k_range = ks_nyq/2;

    % Row-only RF on Nyquist grid (for fair comparison & beamforming)
    v_row1_nyq = zeros(Nt, settings.N_elements);
    v_row2_nyq = zeros(Nt, settings.N_elements);

    % Merged RF on Nyquist grid (THIS replaces the old "overwrite settings" trick)
    v_merge_nyq = zeros(Nt, settings.N_elements);

    %% ----- Build merged interleaved lattice (denser) -----
    % Combined lattice pitch becomes pitch_row/2
    pitch_merge = pitch_row/2;

    % Merged positions (length 2*N_row) sorted left->right
    sp_merge = sort([sp1(:); sp2(:)]).';  % row vector length 2*N_row

    % Mapping indices for placing samples into sp_merge
    [~, idx1] = ismember(sp1(:), sp_merge(:));
    [~, idx2] = ismember(sp2(:), sp_merge(:));

    %% ----- For each time sample: DFT/IDFT onto Nyquist grid -----
    for row = 1:Nt

        % Row-only -> Nyquist
        [Y1, fk] = DFT(v1i(row,:), pitch_row, sp1, k_range);
        [Y2, ~ ] = DFT(v2i(row,:), pitch_row, sp2, k_range);

        v_row1_nyq(row,:) = ks_nyq * IDFT(Y1, fk, sp_nyq, ks_nyq);
        v_row2_nyq(row,:) = ks_nyq * IDFT(Y2, fk, sp_nyq, ks_nyq);

        % Merge channels on denser lattice (no de-alias)
        v_merge = zeros(1, numel(sp_merge));
        v_merge(idx1) = v1i(row,:);
        v_merge(idx2) = v2i(row,:);

        % Merged lattice -> Nyquist
        [Ym, ~] = DFT(v_merge, pitch_merge, sp_merge, k_range);
        v_merge_nyq(row,:) = ks_nyq * IDFT(Ym, fk, sp_nyq, ks_nyq);
    end

    clear v1i v2i Y1 Y2 Ym v_merge;

    %% ----- Group delay compensation -----
    t_offset = round((length(settings.excitation) + ...
                      length(settings.impulse_response) + ...
                      length(settings.impulse_response) - 2)/2) / settings.fs;

    t = (t1 + t2)/2;

    %% ----- Beamforming (all on the SAME Nyquist grid) -----
    bf_row1  = cpwc_beamform_das(v_row1_nyq,  t, t_offset, settings, steer_deg, zero_pad);
    bf_row2  = cpwc_beamform_das(v_row2_nyq,  t, t_offset, settings, steer_deg, zero_pad);
    bf_merge = cpwc_beamform_das(v_merge_nyq, t, t_offset, settings, steer_deg, zero_pad);

    %% ----- Compounding -----
    compound_row1  = compound_row1  + bf_row1;
    compound_row2  = compound_row2  + bf_row2;
    compound_merge = compound_merge + bf_merge;

    clear v_row1_nyq v_row2_nyq v_merge_nyq bf_row1 bf_row2 bf_merge;
end

%% ---------- Image process ----------
img_row1  = img_process(compound_row1',  DR);
img_row2  = img_process(compound_row2',  DR);
img_merge = img_process(compound_merge', DR);

%% ---------- Display (3 panels) ----------
figure('Units','pixels','Position',[100 100 1600 800]);

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile(1)
imagesc(settings.img_x*1e3, settings.img_z*1e3, img_row1);
axis image; colormap(gray(DR)); colorbar;
title('Row 1 only'); xlabel('x [mm]'); ylabel('z [mm]');

nexttile(3)
imagesc(settings.img_x*1e3, settings.img_z*1e3, img_row2);
axis image; colormap(gray(DR)); colorbar;
title('Row 2 only'); xlabel('x [mm]'); ylabel('z [mm]');

nexttile(2, [2 1])
imagesc(settings.img_x*1e3, settings.img_z*1e3, img_merge);
axis image; colormap(gray(DR)); colorbar;
title('Interleaved merge (no de-alias)'); xlabel('x [mm]'); ylabel('z [mm]');

sgtitle(sprintf('CPWC | N=%g | 2-row interleaved, no de-alias', N_row));
saveas(gcf, 'CPWC_row1_row2_merge_no_dealias.fig');




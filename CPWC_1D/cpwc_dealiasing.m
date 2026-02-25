% CPWC_DEALIASING Rigorous CPWC with plotting (script-like behavior)
%
% Usage:
%   cpwc_run();                          % run + plot
%   out = cpwc_run('doPlot',false);      % run only, no plot
%
% This version preserves the original timing logic:
%   tof = delay - t + t_offset
%   pos = round(tof * itp_fs) + zero_pad
clc; clear; close all;
% profile on;
% tic;

settings        = S_DRDP_1D;
compound_deg    = -3:1:3;
point_position  = [0 0 30e-3];
point_amplitude = 1;
zero_pad        = 500;
doPlot          = true;

%% ---------- Field II init ----------
field_init(0);
set_sampling(settings.fs);
set_field('c',  settings.c);
set_field('fs', settings.fs);

%% ---------- CPWC ----------
Nx = numel(settings.img_x);
Nz = numel(settings.img_z);

compound_image = zeros(Nx, Nz);
for k = 1:numel(compound_deg)
    %%  Transducer Setting (fs1)
    steer_deg = compound_deg(k);
    ratio= 3/4;
    [N_ele_1, width_1, pitch_1, sp_1]= transducer_sym_width(settings.N_elements, ratio, settings.kerf, settings.aperture_size);
    
    [v1, t1] = cpwc_one_angle( ...
        settings, steer_deg, point_position, point_amplitude, N_ele_1, width_1, sp_1);
    
    clear ratio;

    %%  Transducer Setting (fs2)
    ratio= 1/2;
    [N_ele_2, width_2, pitch_2, sp_2]= transducer_sym_width(settings.N_elements, ratio, settings.kerf, settings.aperture_size);
    
    [v2, t2] = cpwc_one_angle( ...
        settings, steer_deg, point_position, point_amplitude, N_ele_2, width_2, sp_2);
    
    v2 = v2/(width_2/width_1);
    v1= [v1; zeros(size(v2,1)-size(v1,1), size(v1, 2))];  % Nsamples 補齊（若有需要）
    v1_interp= z_interp(v1, settings.fs, settings.itp_ratio, zero_pad);
    v2_interp= z_interp(v2, settings.fs, settings.itp_ratio, zero_pad);
    v_final= zeros(size(v2_interp, 1), settings.N_elements);
    
    tic;
    %%  De-alias operation
    for row= 1: size(v2_interp, 1)
        
        sp= settings.element_position;    %   space position (Nyquist)
        ks= 1/settings.pitch;             %   space sampling rate (Nyquist)
        k_range= ks/2;                    %   f_range= f_H= f_s/2
    
        %%  DFT of k-space
        %   [X, f]= DFT(sig, Ts, t_axis, f_range);
        [X1, f1]= DFT(v1_interp(row, :), pitch_1, sp_1, k_range);
        [X2, f2]= DFT(v2_interp(row, :), pitch_2, sp_2, k_range);
        % ks/2
        % 1/pitch_2
        %%  IDFT of operation result
        %   de_alias_spec= de_alias_v2(X1, X2, fs, fs1, f)
        R= de_alias_pk(X1, X2, 1/pitch_1, 1/pitch_2, f1);
        x_R= ks*IDFT(R, f2, sp, ks);
        v_final(row, :)= x_R;
    end
    
    toc;

    clear v1_interp v2_interp X1 X2 R x_R;

    %% --- Group-delay compensation (match your original script) ---
    num = (length(settings.excitation) + length(settings.impulse_response) + length(settings.impulse_response) - 2);
    % t_offset = round(num/2) / settings.fs;
    t_offset= round((length(settings.excitation)+length(settings.impulse_response)+length(settings.impulse_response)-2)/2)/settings.fs;
    % tic;
    %%  Beamforming (fast DAS, still rigorous timing)
    t = (t1+t2)/2;
    beamform_image = cpwc_beamform_das(v_final, t, t_offset, settings, steer_deg, zero_pad);
    % toc;
    %%  Compounding image
    compound_image = compound_image+ beamform_image;
    clear v_final reshape_image beamform_image;
end
% profile viewer;

% toc;

% fprintf('compound_image min/max = %.3e / %.3e\n', ...
%     min(compound_image(:)), max(compound_image(:)));
%% ---------- Image processing ----------
DR = settings.dynamic_range;
log_image = img_process(compound_image', DR);

%% ---------- Plot ----------
%   Image Formation
figure;
image(settings.img_x*1e3, settings.img_z*1e3, log_image); % Setting true axis[mm], rotate the log_image
title("De-aliasing before alignment, CPWC");
axis image;
colormap(gray(DR));                                        % Set the gray colormap with dynamic range
cbar= colorbar;                                            % Set the colorbar
title(cbar, 'dB');
xlabel('x[mm]');
ylabel('z[mm]');
saveas(gcf,'CPWC.fig');

%   Beam pattern
figure;
beam_pattern= max(log_image,[], 1);
beam_pattern= beam_pattern-max(beam_pattern);
plot(settings.img_x*1e3, beam_pattern);
title("Beam pattern");
axis tight;
saveas(gcf,'Beam_pattern.fig');

%   Beam pattern (axial)
figure;
beam_pattern_z= max(log_image, [], 2);
beam_pattern_z= beam_pattern_z- max(beam_pattern_z);
plot(settings.img_z*1e3, beam_pattern_z);
title("Beam pattern (Axial)");
axis tight;
saveas(gcf,'Beam_pattern_(Axial).fig');


%% ---------- Output ----------
% out = struct();
% out.settings        = settings;
% out.compound_deg    = compound_deg;
% out.point_position  = point_position;
% out.point_amplitude = point_amplitude;
% out.zero_pad        = zero_pad;
% out.compound_image  = compound_image;
% out.log_image       = log_image;



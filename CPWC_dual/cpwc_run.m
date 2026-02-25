% CPWC_RUN Rigorous CPWC with plotting (script-like behavior)
%
% Usage:
%   cpwc_run();                          % run + plot
%   out = cpwc_run('doPlot',false);      % run only, no plot
%
% This version preserves the original timing logic:
%   tof = delay - t + t_offset
%   pos = round(tof * itp_fs) + zero_pad
clc; clear; close all;
tic;

settings        = S_DRDP_dual;
compound_deg    = -3:1:3;
point_position  = [0 0 30e-3];
point_amplitude = 1;
zero_pad        = 3000;
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
    steer_deg = compound_deg(k);
    % beamform_image 
    
    % If have to change the ratio of element 
    % ratio = 48/64;
    % [N_ele, width, pitch, sp]= transducer_sym_width(settings.N_elements, ratio, settings.kerf, settings.aperture_size);
 
    [v, t] = cpwc_one_angle_dual( ...
        settings, steer_deg, point_position, point_amplitude);%, N_ele, width, sp);  % If have to change the ratio of element 
    v_interp= z_interp(v, settings.fs, settings.itp_ratio, zero_pad);
    %% --- Group-delay compensation (match your original script) ---
    t_offset= round((length(settings.excitation)+length(settings.impulse_response)+length(settings.impulse_response)-2)/2)/settings.fs;
    % tic;
    %%  Beamforming (fast DAS, still rigorous timing)
    beamform_image = cpwc_beamform_das_dual(v_interp, t, t_offset, settings, steer_deg, zero_pad);%, N_ele, sp);  % If have to change the ratio of element 
    % toc;
    %%  Compounding image
    compound_image = compound_image+ beamform_image;
    clear v_final reshape_image beamform_image;
end

toc;

% fprintf('compound_image min/max = %.3e / %.3e\n', ...
%     min(compound_image(:)), max(compound_image(:)));
%% ---------- Image processing ----------
DR = settings.dynamic_range;
log_image = img_process(compound_image', DR);

%% ---------- Plot ----------
if doPlot
    figure;
    imagesc(settings.img_x*1e3, settings.img_z*1e3, log_image);
    axis image; colormap gray; colorbar;
    caxis([0 DR]);
    xlabel('Lateral (mm)');
    ylabel('Depth (mm)');
    title('CPWC (fast, rigorous timing)');
end

%   Beam pattern
figure;
beam_pattern= max(log_image,[], 1);
beam_pattern= beam_pattern-max(beam_pattern);
plot(settings.img_x*1e3, beam_pattern);
title("Beam pattern");
axis tight;
saveas(gcf,'Beam_pattern.fig');

%% ---------- Output ----------
out = struct();
out.settings        = settings;
out.compound_deg    = compound_deg;
out.point_position  = point_position;
out.point_amplitude = point_amplitude;
out.zero_pad        = zero_pad;
out.compound_image  = compound_image;
out.log_image       = log_image;



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

settings        = S_DRDP_dual;
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
    steer_deg = compound_deg(k);
    [v, t] = cpwc_one_angle_dual( ...
        settings, steer_deg, point_position, point_amplitude);
    
    %%  Resample
    %   Nyquist data
    up= 1;
    down= 2;
    ratio= up/down;
    v_rsp= pw_resample(v, up, down);
    [N_ele, width, pitch, sp]= transducer_sym_width(settings.N_elements_x, ratio, settings.kerf_x, settings.aperture_x);
    clear v;

    %   Alias data
    for i= 1: settings.N_elements_y
        %   Odd row (3/4 alias)
        if mod(i, 2)==1
            up= 22;
            down= 64;
            [N_ele_1, width_1, pitch_1, sp_1]= transducer_sym_width(N_ele, up/down, settings.kerf_x, settings.aperture_x);
            v_rsp_fs1(:, :, (i+1)/2)= pw_resample(v_rsp(:, N_ele*(i-1)+1: N_ele*i), up, down);
        
        %   Even row (1/2 alias)
        else
            up= 16;
            down= 64;
            [N_ele_2, width_2, pitch_2, sp_2]= transducer_sym_width(N_ele, up/down, settings.kerf_x, settings.aperture_x);
            v_rsp_fs2(:, :, i/2)= pw_resample(v_rsp(:, N_ele*(i-1)+1: N_ele*i), up, down);
        end
    end
    m= 1; n= 1;
    final_image= zeros(length(settings.img_x), length(settings.img_z));
    for j= 1: settings.N_elements_y-1
        %   Obtain fs1, fs2
        v1= v_rsp_fs1(:, :, m);
        v2= v_rsp_fs2(:, :, n);
        %   Interpolation
        zero_pad= 500;
        v1_interp= z_interp(v1, settings.fs, settings.itp_ratio, zero_pad);
        v2_interp= z_interp(v2, settings.fs, settings.itp_ratio, zero_pad);
    
        %%  De-alias operation
        for row= 1: size(v1_interp, 1)
            % row
            sp= sp;
            ks1= 1/pitch;
            k_range= ks1/2;
        
            %%  DFT of k-space
            %   [X, f]= DFT(sig, fs, t_axis, f_range);
            [X1, f1]= DFT(v1_interp(row, :), pitch_1, sp_1, k_range);
            [X2, f2]= DFT(v2_interp(row, :), pitch_2, sp_2, k_range);
        
            %%  IDFT of operation result
            R= de_alias_pk(X1, X2, 1/pitch_1, 1/pitch_2, f1);
            x_R= ks1*IDFT(R, f1, sp, ks1);
            v_final(row, :)= x_R;
        end
        clear v1_interp v2_interp X1 X2 R x_R;
    
        %%  Beamforming
        t_offset= round((length(settings.excitation)+length(settings.impulse_response)+length(settings.impulse_response)-2)/2)/settings.fs;
        
        %%  Beamforming (fast DAS, still rigorous timing)
        beamform_image = cpwc_beamform_das_dual(v_final, t, t_offset, settings, steer_deg, zero_pad, N_ele, sp);
        final_image= final_image+ beamform_image;
        clear beamform_image;

        %%  Alternative order
        if mod(j, 2)==1
            m= m+1;
        else
            n= n+1;
        end
    end
    compound_image= compound_image+ final_image;
    clear v_rsp_fs1 v_rsp_fs2 final_image
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



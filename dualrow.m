%% simulation_PW1d.m
% Run kN / ks1 / ks2 together and display 3 images

clear; close all; clc;
t_all = tic;

CASE_LIST = {'kN','ks1','ks2'};
DR = 50;

% paper-like view
ov.img_x_border = [-20e-3 20e-3];
ov.img_z_border = [ 20e-3 35e-3];

% phantom
point_position  = [0 0 30e-3];
point_amplitude = 1;

results = struct();

for ic = 1:numel(CASE_LIST)
    CASE = CASE_LIST{ic};

    settings = set_para.get(CASE, ov);

    field_init(0);
    set_field('c', settings.c);
    set_field('fs', settings.fs);

    xe = settings.element_position;

    Th = xdc_linear_array(settings.N_elements, ...
                          settings.width, settings.height, settings.kerf, ...
                          1, 1, [0 0 1]);
    Rh = xdc_linear_array(settings.N_elements, ...
                          settings.width, settings.height, settings.kerf, ...
                          1, 1, [0 0 1]);

    xdc_excitation(Th, settings.excitation);
    xdc_impulse(Th, settings.impulse_response);
    xdc_impulse(Rh, settings.impulse_response);

    tx_delays = -(0 * xe) / settings.c;
    xdc_focus_times(Th, 0, tx_delays');

    [v, t0] = calc_scat_multi(Th, Rh, point_position, point_amplitude);
    v = v - mean(v,1);

    Nt = size(v,1);
    t_offset = round((numel(settings.excitation) + ...
                      2*numel(settings.impulse_response) - 2)/2) ...
                      / settings.fs;

    t_raw = (0:Nt-1)'/settings.fs;
    t_up  = (0:(Nt*settings.itp_ratio-1))'/settings.itp_fs;
    v_up  = interp1(t_raw, v, t_up, 'linear', 0);

    front = 3000; behind = 5000;
    v_pad = [zeros(front, settings.N_elements); ...
             v_up; ...
             zeros(behind, settings.N_elements)];
    idx = (1:size(v_pad,1))';

    [X,Z] = meshgrid(settings.img_x, settings.img_z);
    img = zeros(numel(settings.img_x), numel(settings.img_z));

    for m = 1:settings.N_elements
        t_tx = Z / settings.c;   % plane wave (theta = 0)
        t_rx = sqrt((X-xe(m)).^2 + Z.^2) / settings.c;
        tof = t_tx + t_rx - t0 + t_offset;
        posf = tof * settings.itp_fs + front + 1;
        sig  = interp1(idx, v_pad(:,m), posf, 'linear', 0);
        img  = img + sig.';
    end

    env = abs(hilbert(img.'));
    BdB = 20*log10(env/max(env(:)) + 1e-12);
    disp_dB = max(BdB, -DR) + DR;
    
    results.(CASE).env = env;   % 給 de-aliasing 用
    results.(CASE).disp = disp_dB;
    results.(CASE).x = settings.img_x;
    results.(CASE).z = settings.img_z;

    xdc_free(Th); xdc_free(Rh);
    field_end;
end

% figure('Color','w');
% tiledlayout(3,1,'TileSpacing','compact');

for ic = 1:3
    CASE = CASE_LIST{ic};
    % nexttile;
    figure('Color','w');
    imagesc(results.(CASE).x*1e3, results.(CASE).z*1e3, results.(CASE).disp);
    colormap(gray); caxis([0 DR]); colorbar;
    xlabel('x [mm]'); ylabel('z [mm]');
    title(CASE);
end


fprintf('TOTAL runtime %.2f s\n', toc(t_all));



% ===== De-aliasing (ks1 + ks2) =====
env1 = results.ks1.env;   % [Nz x Nx]
env2 = results.ks2.env;   % [Nz x Nx]
x_im = results.ks1.x;     % meters
z_im = results.ks1.z;     % meters

opts = struct();
opts.use_exclude_center = true;
opts.exclude_center_mm  = 1.0;   % 把主瓣附近 +/- 1mm 排除（你可以調）
opts.pow = 1;

out = dealias_dualpitch(env1, env2, x_im, z_im, opts);

% 顯示 de-aliased score map
figure('Color','w');
imagesc(x_im*1e3, z_im*1e3, out.score);
colormap(gray); colorbar;
xlabel('x [mm]'); ylabel('z [mm]');
title('De-aliased score map = norm(ks1) .* norm(ks2)');
set(gca,'YDir','normal');

% 顯示每個深度的 xhat(z)（紅線）
hold on;
plot(out.xhat_z*1e3, z_im*1e3, 'r.', 'MarkerSize', 6);
hold off;

% 顯示 global peak（綠色叉）
hold on;
plot(out.global_peak.x*1e3, out.global_peak.z*1e3, 'g+', 'MarkerSize', 12, 'LineWidth', 1.5);
hold off;
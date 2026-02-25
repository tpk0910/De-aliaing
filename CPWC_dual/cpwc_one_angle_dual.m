function [v, t] = cpwc_one_angle_dual(settings, steer_deg, point_position, point_amplitude, N_ele_x, width_x, sp_x)
% CPWC_ONE_ANGLE One angle acquisition + DAS beamforming (fast, no pixel loop).
%
% Dependencies (must be on MATLAB path):
%   S_M3_1D.m (settings)
%   pw_steer_delay.m
%   z_interp.m
%   img_process.m
%   Field II: field_init, xdc_linear_array, xdc_excitation, xdc_impulse,
%            xdc_times_focus, calc_scat_multi, xdc_free

%% Default setting
if nargin < 5 || isempty(N_ele_x)
        N_ele_x = settings.N_elements_x;
end

if nargin < 6 || isempty(width_x)
        width_x = settings.width_x;
end

if nargin < 7 || isempty(sp_x)
        sp_x = settings.ele_pos_x;
end

%% --- Transmit aperture ---
enabled = ones(N_ele_x, settings.N_elements_y);
transmit = xdc_2d_array(N_ele_x, settings.N_elements_y, ...
                           width_x     , settings.width_y, ...
                           settings.kerf_x      , settings.kerf_y, ...
                           enabled, 4, 4, [0, 0, 1e4]);
xdc_excitation(transmit, settings.excitation);
xdc_impulse(transmit, settings.impulse_response);

steer_time = pw_steer_delay(steer_deg, sp_x, settings.c);
steer_time_2d= repmat(steer_time, 1, settings.N_elements_y);
xdc_times_focus(transmit, 0, steer_time_2d);

%% --- Receive aperture ---
receive = xdc_2d_array(N_ele_x, settings.N_elements_y, ...
                          width_x     , settings.width_y, ...
                          settings.kerf_x      , settings.kerf_y, ...
                          enabled, 4, 4, [0, 0, 1e4]);
xdc_impulse(receive, settings.impulse_response);

%% --- Channel data from Field II ---
[v, t] = calc_scat_multi(transmit, receive, point_position, point_amplitude);

%% --- Cleanup Field II apertures ---
xdc_free(transmit);
xdc_free(receive);
end

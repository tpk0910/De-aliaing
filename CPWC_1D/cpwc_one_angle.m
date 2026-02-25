function [v, t] = cpwc_one_angle(settings, steer_deg, point_position, point_amplitude, N_ele, width, sp)
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
if nargin < 5 || isempty(N_ele)
        N_ele = settings.N_elements;
end

if nargin < 6 || isempty(width)
        width = settings.width;
end

if nargin < 7 || isempty(sp)
        sp = settings.element_position;
end

%% --- Transmit aperture ---
transmit = xdc_linear_array(N_ele, width, settings.height, ...
    settings.kerf, 4, 4, [0, 0, 1e4]);
xdc_excitation(transmit, settings.excitation);
xdc_impulse(transmit, settings.impulse_response);

steer_time = pw_steer_delay(steer_deg, sp, settings.c);
xdc_times_focus(transmit, 0, steer_time);

%% --- Receive aperture ---
receive = xdc_linear_array(N_ele, width, settings.height, ...
    settings.kerf, 4, 4, [0, 0, 1e4]);
xdc_impulse(receive, settings.impulse_response);

%% --- Channel data from Field II ---
[v, t] = calc_scat_multi(transmit, receive, point_position, point_amplitude);

%% --- Cleanup Field II apertures ---
xdc_free(transmit);
xdc_free(receive);
end

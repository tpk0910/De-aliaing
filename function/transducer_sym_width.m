%   INTRO: Transducer setting, creating WIDE element
%   PARAM: 
%   (1) pitch (change to pitch)
%   (2) ratio= lower rate, compared to non-aliasing fs
%   (3) kerf
%   (4) aperture_size
% (4.4000e-04, 3/4, 1.0000e-05, 0.0141)
function [N_ele_new, width_new, pitch_new, element_position]= transducer_sym_width(N_ele, ratio, kerf, aperture_size)
    N_ele_new = round((N_ele * ratio) / 2)*2;
    pitch_new = aperture_size*2 / N_ele_new;
    width_new = pitch_new - kerf;
    element_position = (-(aperture_size - (pitch_new/2))) : pitch_new : (aperture_size - (pitch_new/2));
    % fprintf('pitch_new = %d\n', pitch_new);
    % fprintf('element_position = %d\n', element_position);
    % fprintf('width_new = %d\n', width_new);
    % fprintf('N_ele_new = %d\n', N_ele_new);
end

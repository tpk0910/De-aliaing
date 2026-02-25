%   INTRO: Steer plane wave transmit time delay
%   PARAM: 
%   (1) deg= steering angle, degree
%   (2) ele_pos= position of each element
%   (3) c= sound_velocity

function [steer_time_delay]= pw_steer_delay(deg, ele_pos, c)
    %   Focus settings
    theta= deg/180*pi;                       %   Steering Angle (rad)
    dist= ele_pos.*sin(theta);           %   Delay Distance for Each Element
    steer_time_delay= dist/c;                      %   Delay Time for Each Element
end
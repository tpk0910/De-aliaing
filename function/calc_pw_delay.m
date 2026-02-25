%   INTRO: Calculate pixel-wise two way delay of plane wave transmission
%   PARAM: 
%   (1) N_ele= number of elements
%   (2) element_position= position of each element
%   (3) img_x= user decided image grid (x)
%   (4) img_z= user decided image grid (z)
%   (5) c= sound velocity
%   (6) theta= plane wave steering angle, degree

function [delay] = calc_pw_delay(N_ele, ele_pos, img_x, img_z, c, deg)
    %   Convert deg to radian
    theta= deg/180*pi;    
    %   Pixel grid (Store delay value of each pixel)
    delay= zeros(N_ele, length(img_x), length(img_z));
    for i= 1: N_ele
        for j= 1:length(img_x)
            tx_delay= img_x(j)*sin(theta);
            tz_delay= img_z.*cos(theta);
            rx_delay= (ele_pos(i)- img_x(j))^2;
            rz_delay= (img_z).^2;
            delay(i,j,:) = tx_delay+ tz_delay+ sqrt(rx_delay+rz_delay);
        end
    end
    delay = delay/c;
end
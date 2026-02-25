%%  INTRO
%   Resample (interpolation+ summation)

function [v_rsp]= pw_resample(v, up, down)
    x_ori= 1: size(v, 2);
    x_itp= linspace(1, max(x_ori), up*length(x_ori));
    for row= 1: size(v, 1)
        v_itp(row, :)= interp1(x_ori, v(row, :), x_itp, "linear");
    end
    
    %%  Down-sample
    % v_rsp= v_itp(:, 1: down: end);

    %%  Sum
    v_rsp= reshape(v_itp, size(v_itp, 1), down, size(v_itp, 2)/down);
    v_rsp= squeeze(sum(v_rsp, 2))/down;
end
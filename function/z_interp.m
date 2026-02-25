%   INTRO: Interpolation and zero-padding for channel data
%   *< In US imaging system, RF data should be processed to form an image.  
%      (1) Demodulation or (2) Interpolation + Hilbert transform, must 
%      choose one to execute to channel data >*
%   PARAM: 
%   (1) v= channel data
%   (2) fs= sampling frequency
%   (3) ratio= interpolation rate
%   (4) zero_pad= padding number


function v_interp= z_interp(v, fs, ratio, zero_pad)

   t_axis= (0: size(v, 1)-1)/fs;
   t_interp= (0: size(v,1)*ratio-1)/(fs*ratio);
   v_interp= interp1(t_axis, v, t_interp, "linear");
   v_interp(isinf(v_interp))= 0;
   v_interp(isnan(v_interp))= 0;
   
   %%   Check for FDS or RF data
   if (size(v, 3)==1)
       v_interp= [zeros(zero_pad, size(v, 2)); v_interp];
       v_interp= [v_interp; zeros(zero_pad, size(v, 2))];
   else
       v_interp= [zeros(zero_pad, size(v, 2), size(v, 3)); v_interp];
       v_interp= [v_interp; zeros(zero_pad, size(v, 2), size(v, 3))];
   end
end
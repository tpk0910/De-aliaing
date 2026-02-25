%   INTRO: Form an image, and calculate the energy (if needed)
%   PARAM: 
%   (1) bf_image= beamformed_image
%   (2) DR= dynamic_range
%   VARARGIN: 
%   (3) varargin{1}= ROI= region of interest for energy calculation

function [log_image, energy]= img_process(bf_image, DR, varargin)
    
    %   Hilbert Transform
    envelope= abs(hilbert(bf_image));
    envelope(isnan(envelope)|isinf(envelope))= 0;
    
    if nargin >= 3 && nargout >= 2
        %   Energy Calculation
        ROI= logical(varargin{1});
        dx= varargin{2};
        dz= varargin{3};
        pixel= dx*dz;

        power= envelope.^2;
        energy = sum(power(ROI), 'all')*pixel;
    else 
        energy= [];
    end

    %   Log Compression
    %   (Prof, after demodulation)
    nor_image= envelope./max(envelope(:));
    log_image= 20*log10(nor_image)+ DR;
    log_image(log_image<0)= 0;
    
    % %   (Senior, after demodulation)
    % log_image= 20*log10(envelope);
    % log_image= log_image-max(log_image(:))+DR;
    % log_image(log_image<0)= 0;

end
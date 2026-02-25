%%  Intro:
%   Set file for US simulation image
%   Standard pitch
%   Static width

classdef S_DRDP_dual
    properties(Constant)
        
        %% Define the transducer     
        f0 = 3.5e6;                                     % Center frequency [Hz]       
        fs = 14e6;                                      % Sampling frequency [Hz] 
        c = 1540;                                       % Speed of sound [m/s] 
        lambda = S_DRDP_dual.c/S_DRDP_dual.f0;                            % Wavelength [m] 
        N_elements_x = 128;                              % Number of elements in x 
        N_elements_y = 2;                             % Number of elements in y
        pitch_x= S_DRDP_dual.lambda/2*(128/S_DRDP_dual.N_elements_x);                             % Standard pitch size [m]
        pitch_y= S_DRDP_dual.lambda*1;                           % y-direction pitch for 1.5D array
        kerf_x= 0.01e-3;
        kerf_y= 0.01e-3;
        width_x= S_DRDP_dual.pitch_x-S_DRDP_dual.kerf_x;
        width_y= S_DRDP_dual.pitch_y-S_DRDP_dual.kerf_y;
        aperture_x = S_DRDP_dual.pitch_x*S_DRDP_dual.N_elements_x/2;      % Aperture size [m] 
        aperture_y = S_DRDP_dual.pitch_y*S_DRDP_dual.N_elements_y/2;      % Aperture size [m] 
        ele_pos_x= (-(S_DRDP_dual.N_elements_x/2)+0.5)*S_DRDP_dual.pitch_x:...
                                              S_DRDP_dual.pitch_x:...
                    ((S_DRDP_dual.N_elements_x/2)-0.5)*S_DRDP_dual.pitch_x;
        ele_pos_y= (-(S_DRDP_dual.N_elements_y/2)+0.5)*S_DRDP_dual.pitch_y:...
                                              S_DRDP_dual.pitch_y:...
                    ((S_DRDP_dual.N_elements_y/2)-0.5)*S_DRDP_dual.pitch_y; 
        
        %%  Conventional B-mode
        focal_depth= 30e-3;
                
        %% Define the impulse response and excitation pulse of the transducer 
        cycles= 2;
        t= 0: 1/S_DRDP_dual.fs : S_DRDP_dual.cycles/S_DRDP_dual.f0;

        % Using sine wave as excitation function
        excitation= sin(2*pi*S_DRDP_dual.f0*S_DRDP_dual.t);        
        impulse_response= S_DRDP_dual.excitation.*hanning(max(size(S_DRDP_dual.excitation)))'; 
        pulse_length= S_DRDP_dual.cycles*S_DRDP_dual.lambda;
       
        %%  Self-define the image grid (Spatial domain)
        dx= S_DRDP_dual.lambda/4;
        dy= S_DRDP_dual.lambda/4;
        dz= S_DRDP_dual.c/S_DRDP_dual.fs/2;

        img_start= 20; 
        img_end= 40;

        % In order to see the grating lobes, img_x_border > aperture_size 
        img_x_border= [-S_DRDP_dual.aperture_x*1.5 S_DRDP_dual.aperture_x*1.5];      % Lateral border of image grid
        img_y_border= [-S_DRDP_dual.aperture_y*1.5 S_DRDP_dual.aperture_y*1.5];  % Elevational border of image grid
        img_z_border= [S_DRDP_dual.img_start/1000  S_DRDP_dual.img_end/1000];        % Axial border of image grid
        
        img_x= S_DRDP_dual.img_x_border(1):S_DRDP_dual.dx:S_DRDP_dual.img_x_border(2);     % Lateral axis of image
        img_y= S_DRDP_dual.img_y_border(1):S_DRDP_dual.dy:S_DRDP_dual.img_y_border(2);     % Elevational axis of image
        img_z= S_DRDP_dual.img_z_border(1):S_DRDP_dual.dz:S_DRDP_dual.img_z_border(2);     % Axial axis of image

        num_x_pixel= length(S_DRDP_dual.img_x);      % Number of lateral pixel 
        num_y_pixel= length(S_DRDP_dual.img_y);      % Number of elevational pixel
        num_z_pixel= length(S_DRDP_dual.img_z);      % Number of axial pixel

        %% Define the properties of image
        dynamic_range= 50;
        
        %% Define the K-space (Spatial frequency domain)
        % Derived from self-define image grid, 1/(FOV_x)= dkx
        k0= S_DRDP_dual.f0/S_DRDP_dual.c;

        dkx = 1/(S_DRDP_dual.num_x_pixel*S_DRDP_dual.dx);
        kx = [-(0.5/S_DRDP_dual.dx):S_DRDP_dual.dkx:(0.5/S_DRDP_dual.dx-S_DRDP_dual.dkx)]./S_DRDP_dual.k0;

        dkz = 1/(S_DRDP_dual.num_z_pixel*S_DRDP_dual.dz);
        kz = [flip(-(0.5/S_DRDP_dual.dz):S_DRDP_dual.dkz:(0.5/S_DRDP_dual.dz))]./S_DRDP_dual.k0;

        %% Define the baseband
        order = 3;
        fd = S_DRDP_dual.f0;
        fcut = 0.5*S_DRDP_dual.fd;
        
        % Define the channel data                   
        itp_ratio = 8;
        itp_fs = S_DRDP_dual.itp_ratio*S_DRDP_dual.fs;   

    end   
end
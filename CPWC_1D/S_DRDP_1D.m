%%  Intro:
%   Set file for US simulation image
%   Standard pitch
%   Static width

classdef S_DRDP_1D
    properties(Constant)
        
        %% Define the transducer     
        f0 = 3.5e6;                                     % Center frequency [Hz]       
        fs = 14e6;                                      % Sampling frequency [Hz] 
        c = 1540;                                       % Speed of sound [m/s] 
        lambda = S_DRDP_1D.c/S_DRDP_1D.f0;                              % Wavelength [m] 
        N_elements = 64;                               % Number of elements 
        pitch= S_DRDP_1D.lambda/2*(128/S_DRDP_1D.N_elements);           % Standard pitch size [m]
        kerf= 0.01e-3;
        width= S_DRDP_1D.pitch-S_DRDP_1D.kerf;
        height = 12e-3;                                 % Height of element [m] 
        aperture_size = S_DRDP_1D.pitch*S_DRDP_1D.N_elements/2;         % Aperture size [m] (only positive half because of N_elements/2)
        element_position= (-(S_DRDP_1D.N_elements/2)+0.5)*S_DRDP_1D.pitch:...
                                                  S_DRDP_1D.pitch:...
                           ((S_DRDP_1D.N_elements/2)-0.5)*S_DRDP_1D.pitch;  
        %%  Conventional B-mode
        focal_depth= 30e-3;
                
        %% Define the impulse response and excitation pulse of the transducer 
        cycles= 2;
        t= 0: 1/S_DRDP_1D.fs : S_DRDP_1D.cycles/S_DRDP_1D.f0;

        % Using sine wave as excitation function
        excitation= sin(2*pi*S_DRDP_1D.f0*S_DRDP_1D.t);        
        impulse_response= S_DRDP_1D.excitation.*hanning(max(size(S_DRDP_1D.excitation)))'; 
        pulse_length= S_DRDP_1D.cycles*S_DRDP_1D.lambda;
       
        %%  Self-define the image grid (Spatial domain)
        dx= S_DRDP_1D.lambda/4;
        dz= S_DRDP_1D.c/S_DRDP_1D.fs/2;

        img_start= 20; 
        img_end= 40;

        % In order to see the grating lobes, img_x_border > aperture_size 
        img_x_border= [-S_DRDP_1D.aperture_size*1.5 S_DRDP_1D.aperture_size*1.5];       % Lateral border of image grid
        img_z_border= [S_DRDP_1D.img_start/1000 S_DRDP_1D.img_end/1000];        % Axial border of image grid
        
        img_x= S_DRDP_1D.img_x_border(1):S_DRDP_1D.dx:S_DRDP_1D.img_x_border(2);     % Lateral axis of image
        img_z= S_DRDP_1D.img_z_border(1):S_DRDP_1D.dz:S_DRDP_1D.img_z_border(2);     % Axial axis of image

        num_x_pixel= length(S_DRDP_1D.img_x);      % Number of lateral pixel 
        num_z_pixel= length(S_DRDP_1D.img_z);      % Number of axial pixel

        %% Define the properties of image
        dynamic_range= 50;
        
        %% Define the K-space (Spatial frequency domain)
        % Derived from self-define image grid, 1/(FOV_x)= dkx
        k0= S_DRDP_1D.f0/S_DRDP_1D.c;

        dkx = 1/(S_DRDP_1D.num_x_pixel*S_DRDP_1D.dx);
        kx = [-(0.5/S_DRDP_1D.dx):S_DRDP_1D.dkx:(0.5/S_DRDP_1D.dx-S_DRDP_1D.dkx)]./S_DRDP_1D.k0;

        dkz = 1/(S_DRDP_1D.num_z_pixel*S_DRDP_1D.dz);
        kz = [flip(-(0.5/S_DRDP_1D.dz):S_DRDP_1D.dkz:(0.5/S_DRDP_1D.dz-S_DRDP_1D.dkz))]./S_DRDP_1D.k0;

        %% Define the baseband
        order = 3;
        fd = S_DRDP_1D.f0;
        fcut = 0.5*S_DRDP_1D.fd;
        
        % Define the channel data                   
        itp_ratio = 8;
        itp_fs = S_DRDP_1D.itp_ratio*S_DRDP_1D.fs;   

    end   
end
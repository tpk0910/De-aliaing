classdef set_para
    properties(Constant)

        %% ===================== SYSTEM =====================
        f0 = 3.5e6;                         % [Hz]
        fs = 32 * set_para.f0;              % [Hz]
        c  = 1540;                          % [m/s]
        lambda = set_para.c / set_para.f0;  % [m]

        %% ===================== ELEMENT GEOMETRY =====================
        height = 8.5e-3;                    % [m]
        kerf   = 0.01e-3;                   % [m]

        %% ===================== PULSE =====================
        cycles = 2;

        %% ===================== IMAGE GRID DEFAULT =====================
        dx = set_para.lambda/4;
        dz = set_para.c/set_para.fs/2;

        %% ===================== INTERPOLATION =====================
        itp_ratio = 8;
        itp_fs = set_para.itp_ratio * set_para.fs;

        %% ===================== ARRAY CASES (paper table) =====================
        % kN : d = 1*lambda,  Nele=64
        kN_Nele  = 64;
        kN_pitch = 1.0 * set_para.lambda;
        kN_width = set_para.kN_pitch - set_para.kerf;

        % ks1: d = 4/3*lambda, Nele=48
        ks1_Nele  = 48;
        ks1_pitch = (4/3) * set_para.lambda;
        ks1_width = set_para.ks1_pitch - set_para.kerf;

        % ks2: d = 2*lambda,  Nele=32
        ks2_Nele  = 32;
        ks2_pitch = 2.0 * set_para.lambda;
        ks2_width = set_para.ks2_pitch - set_para.kerf;
    end

    methods(Static)
        function settings = get(case_name, overrides)
            if nargin < 2
                overrides = struct();
            end

            %% ---- system ----
            settings.f0 = set_para.f0;
            settings.fs = set_para.fs;
            settings.c  = set_para.c;
            settings.lambda = set_para.lambda;

            settings.height = set_para.height;
            settings.kerf   = set_para.kerf;

            %% ---- pulse ----
            settings.cycles = set_para.cycles;
            settings.t = 0:1/settings.fs:settings.cycles/settings.f0;
            settings.excitation = sin(2*pi*settings.f0*settings.t);
            settings.impulse_response = settings.excitation .* ...
                                        hanning(numel(settings.excitation))';
            settings.pulse_length = settings.cycles * settings.lambda;

            %% ---- interpolation ----
            settings.itp_ratio = set_para.itp_ratio;
            if isfield(overrides,'itp_ratio')
                settings.itp_ratio = overrides.itp_ratio;
            end
            settings.itp_fs = settings.itp_ratio * settings.fs;

            %% ---- choose case ----
            switch lower(case_name)
                case 'kn'
                    settings.case_name  = 'kN';
                    settings.N_elements = set_para.kN_Nele;
                    settings.pitch      = set_para.kN_pitch;
                    settings.width      = set_para.kN_width;
                case 'ks1'
                    settings.case_name  = 'ks1';
                    settings.N_elements = set_para.ks1_Nele;
                    settings.pitch      = set_para.ks1_pitch;
                    settings.width      = set_para.ks1_width;
                case 'ks2'
                    settings.case_name  = 'ks2';
                    settings.N_elements = set_para.ks2_Nele;
                    settings.pitch      = set_para.ks2_pitch;
                    settings.width      = set_para.ks2_width;
                otherwise
                    error('Unknown case_name. Use kN / ks1 / ks2');
            end

            %% ---- derived geometry ----
            settings.aperture_size = settings.pitch * settings.N_elements;
            settings.element_position = (-(settings.N_elements/2)+0.5)*settings.pitch : ...
                                         settings.pitch : ...
                                        ((settings.N_elements/2)-0.5)*settings.pitch;
            settings.element_position = settings.element_position(:); % column

            %% ---- image grid ----
            settings.dx = set_para.dx;
            settings.dz = set_para.dz;

            % default paper-ish borders if not overridden
            settings.img_x_border = [-20e-3 20e-3];
            settings.img_z_border = [ 20e-3 35e-3];

            if isfield(overrides,'img_x_border'), settings.img_x_border = overrides.img_x_border; end
            if isfield(overrides,'img_z_border'), settings.img_z_border = overrides.img_z_border; end
            if isfield(overrides,'dx'), settings.dx = overrides.dx; end
            if isfield(overrides,'dz'), settings.dz = overrides.dz; end

            settings.img_x = settings.img_x_border(1):settings.dx:settings.img_x_border(2);
            settings.img_z = settings.img_z_border(1):settings.dz:settings.img_z_border(2);
        end
    end
end

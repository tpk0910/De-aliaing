%   INTRO: Inverse discrete Fourier transform (IDFT)
%   PARAM: 
%   (1) X= DFT spectrum
%   (2) f= frequency axis
%   (3) t= time axis
%   (4) fs= sampling frequency

function x= IDFT(X, f, t, fs)

    M= length(f);                                   %   頻域總點數
    m= -(M-mod(M,2))/2 : (M-mod(M,2))/2;            %   頻域整數序列（也須對稱）
    N= length(t);                                   %   時域總點數
    n= t*fs;                                        %   時域整數序列（也須對稱）
    x= X*exp(j*2*pi*m'/M*n);
    x= x/M;

end
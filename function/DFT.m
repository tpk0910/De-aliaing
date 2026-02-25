%   INTRO: Discrete Fourier transform (DFT)
%   PARAM: 
%   (1) x_func= input signal
%   (2) Ts= sampling interval
%   (3) t= time axis
%   (4) f_range= desired DFT spectrum range

function [X, f]= DFT(x_func, Ts, t, f_range)

    %%  時域設定
    dt= Ts;                                     %   時間取樣間隔
    N= length(t);                               %   時域總點數
    n= t/Ts;                                    %   時域整數序列
    L= N*dt;                                %   時間長度

    %%  頻域設定
    %   當取樣率變低=取樣間隔變大=點數變少，但時間長度L=N*Ts是固定的
    %   確保不同取樣率的訊號，能在頻域上有相同的解析度

    df= 1/L;                                    %   頻率取樣間隔
    fr= [0: df: f_range];                       %   右側頻率軸（因需以0為中心往左右對稱延伸，
    fl= [-max(fr): df: 0-df];                   %   左側頻率軸　因此左右軸各自設定好再合併）
    f= [fl, fr];                                %   頻率軸
    M= length(f);                               %   頻域總點數
    m= -(M-1)/2 : (M-1)/2;                      %   頻域整數序列（也須對稱）

    %%  訊號與DFT操作
    x= x_func;                                  %   主函數
    X= (x*exp(-j*2*pi*n'/N*m));                 %   DFT操作
    X= X*Ts;
end

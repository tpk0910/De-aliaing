%   INTRO: De-aliasing process
%   PARAM: 
%   (1) X1= DFT spectrum (Undersampled but higher fs)
%   (2) X2= DFT spectrum (Undersampled, 1/2 ideal fs)
%   (3) fs= Nyquist (Ideal) sampling frequency
%   (4) fs1= sampling frequency of X1 (undersampled by higher fs)
%   (5) f= frequency axis

function de_alias_spec= de_alias_v2(X1, X2, fs, fs1, f)
    
    interval= fs1-fs/2;                 % 區段長度
    % fs/2/interval
    iteration= round(fs/2/interval)-1;  % 迴圈次數
    % iteration
    for i= 0: iteration-1
        j= i+1;                                 % i= 左側索引, j= 右側索引
        oper_sector= zeros(1, length(X1));      % 正在執行的區段
        oper_sector((f>= i*interval)& (f<= j*interval))=1;
    
        oper_X1= X1.*oper_sector;           % fs1 頻譜正在執行區域
        oper_X2= X2.*oper_sector;           % fs2 頻譜正在執行區域
    
        subtraction= oper_X2-oper_X1;       % 相減
        df= f(2)-f(1);                      % 頻率間隔
        shift_point= round(interval/df);    % 預計移動點數（長度／間隔）
        shifting= circshift(subtraction, shift_point);  % 實際移動
        flipping= conj(flip(shifting, 2));  % 另一側的混疊區
        symmetry= shifting+ flipping;       % 兩側的混疊區
        
        X1= X1-symmetry;                    % 部分清除的 fs1 頻譜
    end
    de_alias_spec= X1;                      % 完整清除的 fs1 頻譜 
end

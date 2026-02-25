%   INTRO: De-aliasing process
%   PARAM: 
%   (1) X1= DFT spectrum (Undersampled but higher fs)
%   (2) X2= DFT spectrum (Undersampled, 1/2 ideal fs)
%   (3) fs1= sampling frequency of X1 (undersampled by higher fs)
%   (4) fs2= sampling frequency of X2 (undersampled by lower fs)
%   (5) f= frequency axis
% (4/3, 1/2)
% interval
function de_alias_spec= de_alias_pk(X1, X2, fs1, fs2, f)
    
    N = numel(f);                       % f總共有幾個點
    interval = fs1 - fs2;                  % 區段長度(Δf)
    iteration = round(fs2/interval)-1;   % 迴圈次數
    iteration = max(iteration, 1);      % 確保迴圈至少有一次
    df = f(2) - f(1); % 每個頻率點的間距
    shift_point = round(interval / df); % 把區間從頻率變成index(一次位移幾個點)
    
    mirror_idx = N:-1:1; % 對應到負頻率的索引

    % (4/3, 1/2)
    % interval=4/3-1/2=8-3/6=5/6
    % iteration
    for i = 0 : iteration-1               
        oper_section = (f >= (i * interval)) & (f < ((i+1) * interval)); % i= 左側索引, i+1= 右側索引
        
        if ~any(oper_section), continue; end % 如果不是在f_range裡的話就跳過不處理

        idx = find(oper_section); % 正在執行的區段 (flag)
        % oper_section可能會像是
        % idx = find(oper_section)是找裡面有執行的區段拿出來，不然直接對oper_section會出錯
        % ex: oper_section = [0 0 1 1 0 0]
        %     find(oper_section) = [3 4]

        sub = X2(idx) - X1(idx); % 相減
        
        % 相當於circshift但只有對要處理的區段做shift, circshift是對整個向量
        idx_shift = idx + shift_point; 
        idx_shift = mod(idx_shift-1, N) + 1;
    
        idx_m = mirror_idx(idx_shift); % 找另外一側的
         
        all_idx = [idx_shift(:); idx_m(:)];      % 合併兩側要處理的index
        all_sub = [-sub(:);      -conj(sub(:))]; % 合併兩側要處理的混疊區

        delta = accumarray(all_idx, all_sub, [N 1], @sum, 0); % 把指定index要減掉的地方對好

        X1 = X1 + reshape(delta, size(X1)); % 部分清除的 fs1 頻譜
    end
    de_alias_spec = X1; % 完整清除的 fs1 頻譜
end

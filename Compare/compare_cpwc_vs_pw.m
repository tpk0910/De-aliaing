function T = compare_cpwc_vs_pw(cfg)
%COMPARE_CPWC_VS_PW
% Robust compare between your cpwc_dealiasing.m and senior's PW_dealiasing.m
% even though they are in different folders and begin with:
%   clc; clear; close all;
%
% What it does:
% - Runs each script in an isolated function workspace (so "clear" won't kill cfg)
% - Captures the 3 figures produced by each script
%   Fig1: 2D image (B-mode) -> compute FWHM/PSLR/ISLR on 2D PSF proxy
%   Fig2: 1D lateral beam pattern (plot) -> compute FWHM/PSLR/ISLR on 1D
%   Fig3: 1D axial beam pattern (plot) -> compute FWHM/PSLR/ISLR on 1D
% - Generates overlay/diff/profiles and a summary table.
%
% Usage:
% cfg = struct;
% cfg.cpwc_dir   = "D:\...\cpwc_folder";
% cfg.pw_dir     = "D:\...\pw_folder";
% cfg.cpwc_entry = "cpwc_dealiasing";   % script name without .m
% cfg.pw_entry   = "PW_dealiasing";
% cfg.save_dir   = "D:\...\compare_out"; % "" to disable saving
% cfg.showFigures = true;
% T = compare_cpwc_vs_pw(cfg);

% ---------------- defaults ----------------
if ~isfield(cfg,'cpwc_entry'), cfg.cpwc_entry = 'cpwc_dealiasing'; end
if ~isfield(cfg,'pw_entry'),   cfg.pw_entry   = 'PW_dealiasing';   end
if ~isfield(cfg,'expected_figs'), cfg.expected_figs = 3; end
if ~isfield(cfg,'alpha'), cfg.alpha = 0.5; end
if ~isfield(cfg,'save_dir'), cfg.save_dir = ''; end
if ~isfield(cfg,'showFigures'), cfg.showFigures = true; end

if ~isfolder(cfg.cpwc_dir), error("cpwc_dir not found: %s", cfg.cpwc_dir); end
if ~isfolder(cfg.pw_dir),   error("pw_dir not found: %s", cfg.pw_dir); end

if ~isempty(cfg.save_dir) && ~isfolder(cfg.save_dir)
    mkdir(cfg.save_dir);
end

% ---------------- run + capture ----------------
Me = run_and_capture_three(cfg.cpwc_dir, cfg.cpwc_entry, cfg.expected_figs);
MingChi = run_and_capture_three(cfg.pw_dir,   cfg.pw_entry,   cfg.expected_figs);

% ---------------- compare pairwise ----------------
rows = [];
pairNames = "Fig" + string(1:cfg.expected_figs);

for k = 1:cfg.expected_figs
    a = Me{k}; b = MingChi{k};

    if a.kind == "image2d" && b.kind == "image2d"
        out = compare_image2d(a, b, k, cfg.alpha, cfg.showFigures);
        rows = [rows; {pairNames(k), "2D", out.corr, out.ssim, out.mse, out.psnr_db, out.mad, out.absdiff_p95, ...
            out.FWHM_lat_A, out.FWHM_lat_B, out.PSLR_dB_A, out.PSLR_dB_B, out.ISLR_dB_A, out.ISLR_dB_B}];
        if ~isempty(cfg.save_dir)
            save_all_figs(out, cfg.save_dir, pairNames(k));
        end

    elseif a.kind == "line1d" && b.kind == "line1d"
        out = compare_line1d(a, b, k, cfg.showFigures);
        rows = [rows; {pairNames(k), "1D", out.corr, NaN, out.mse, out.psnr_db, out.mad, out.absdiff_p95, ...
            out.FWHM_A, out.FWHM_B, out.PSLR_dB_A, out.PSLR_dB_B, out.ISLR_dB_A, out.ISLR_dB_B}];
        if ~isempty(cfg.save_dir)
            save_all_figs(out, cfg.save_dir, pairNames(k));
        end

    else
        % mixed types -> fallback to screenshot compare (still works)
        warning("Fig %d is mixed type (%s vs %s). Falling back to screenshot image compare.", ...
            k, a.kind, b.kind);
        aa = ensure_image2d(a);
        bb = ensure_image2d(b);
        out = compare_image2d(aa, bb, k, cfg.alpha, cfg.showFigures);
        rows = [rows; {pairNames(k), "Mixed->2D", out.corr, out.ssim, out.mse, out.psnr_db, out.mad, out.absdiff_p95, ...
            out.FWHM_lat_A, out.FWHM_lat_B, out.PSLR_dB_A, out.PSLR_dB_B, out.ISLR_dB_A, out.ISLR_dB_B}];
        if ~isempty(cfg.save_dir)
            save_all_figs(out, cfg.save_dir, pairNames(k));
        end
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'Pair','Kind','Corr','SSIM','MSE','PSNR_dB','MAD','AbsDiff_P95', ...
     'FWHM_A','FWHM_B','PSLR_dB_A','PSLR_dB_B','ISLR_dB_A','ISLR_dB_B'});

disp("=== Summary table ===");
disp(T);

if ~isempty(cfg.save_dir)
    writetable(T, fullfile(cfg.save_dir, "compare_summary.csv"));
end

end

% ======================================================================
%                           RUN + CAPTURE
% ======================================================================

function contents = run_and_capture_three(folder, entry, expectedN_in)
% Bullet-proof runner against scripts that call: clc; clear; close all;
% We execute the script in BASE workspace so its "clear" doesn't wipe our local vars.

folder = char(folder);
entry  = char(entry);
Nfig   = double(expectedN_in);   % local constant (function workspace)

% save current env
p0 = path;
c0 = pwd;

% temp dir to avoid writing into project folder
tmp = fullfile(tempdir, ['compare_run_tmp_' char(java.util.UUID.randomUUID)]);
mkdir(tmp);

% --- run in BASE workspace ---
evalin('base', 'close all;');  % clean figures in base
evalin('base', sprintf('addpath(''%s'');', folder));
evalin('base', sprintf('cd(''%s'');', tmp));

% IMPORTANT: run the script in base
scriptPath = fullfile(folder, [entry '.m']);
evalin('base', sprintf('run(''%s'');', scriptPath));

% --- NOW fetch figures from ROOT (not base vars) ---
figs = findall(0,'Type','figure');
figs = sort_figs_by_number(figs);

if numel(figs) < Nfig
    % restore
    path(p0); cd(c0);
    error('Script %s produced %d figures (expected %d).', entry, numel(figs), Nfig);
end

figs = figs(1:Nfig);
contents = cell(1,Nfig);
for k = 1:Nfig
    contents{k} = extract_primary_content(figs(k));
end

% restore
path(p0);
cd(c0);

end


function figs = sort_figs_by_number(figs)
nums = arrayfun(@(f) f.Number, figs);
[~,idx] = sort(nums);
figs = figs(idx);
end

% ======================================================================
%                    CONTENT EXTRACTION (image or line)
% ======================================================================

function S = extract_primary_content(fig)
% Pick the "main axes" (largest area, not colorbar), then extract:
% - image2d: matrix + axes scales
% - line1d: x,y data
% - fallback: screenshot image2d

ax = findall(fig, 'Type','axes');
ax = ax(~strcmpi(get(ax,'Tag'),'Colorbar')); % remove colorbar axes if tagged
ax = ax(~strcmpi(get(ax,'Type'),'legend'));  % defensive

if isempty(ax)
    % fallback screenshot
    fr = getframe(fig);
    S.kind = "image2d";
    S.I = im2double(fr.cdata);
    S.x = 1:size(S.I,2);
    S.z = 1:size(S.I,1);
    S.name = string(fig.Number);
    return;
end

% choose axes with max area
areas = zeros(size(ax));
for i=1:numel(ax)
    pos = ax(i).Position;
    areas(i) = pos(3)*pos(4);
end
[~,ii] = max(areas);
ax0 = ax(ii);

% 1) image object?
hImg = findall(ax0,'Type','image');
if ~isempty(hImg)
    I = hImg(1).CData;
    S.kind = "image2d";
    S.I = im2double(I);
    % x/y scale if available
    try
        xd = hImg(1).XData; yd = hImg(1).YData;
        if numel(xd)==2, S.x = linspace(xd(1),xd(2), size(S.I,2)); else, S.x = xd; end
        if numel(yd)==2, S.z = linspace(yd(1),yd(2), size(S.I,1)); else, S.z = yd; end
    catch
        S.x = 1:size(S.I,2); S.z = 1:size(S.I,1);
    end
    S.name = string(fig.Number);
    return;
end

% 2) surface object? (pcolor/surf)
hSurf = findall(ax0,'Type','surface');
if ~isempty(hSurf)
    I = hSurf(1).CData;
    S.kind = "image2d";
    S.I = im2double(I);
    S.x = 1:size(S.I,2); S.z = 1:size(S.I,1);
    S.name = string(fig.Number);
    return;
end

% 3) line plot?
hLine = findall(ax0,'Type','line');
if ~isempty(hLine)
    x = hLine(1).XData(:);
    y = hLine(1).YData(:);
    S.kind = "line1d";
    S.x = double(x);
    S.y = double(y);
    S.name = string(fig.Number);
    return;
end

% 4) fallback screenshot of axes
fr = getframe(ax0);
S.kind = "image2d";
S.I = im2double(fr.cdata); % RGB
S.x = 1:size(S.I,2);
S.z = 1:size(S.I,1);
S.name = string(fig.Number);
end

function S2 = ensure_image2d(S)
if S.kind == "image2d"
    S2 = S; return;
end
% line -> render to image by plotting into an offscreen fig then getframe
f = figure('Visible','off');
plot(S.x, S.y); grid on;
fr = getframe(gca);
close(f);
S2.kind = "image2d";
S2.I = im2double(fr.cdata);
S2.x = 1:size(S2.I,2);
S2.z = 1:size(S2.I,1);
S2.name = S.name;
end

% ======================================================================
%                          2D IMAGE COMPARE
% ======================================================================

function out = compare_image2d(A, B, figIndex, alpha, showFigures)
% grayscale normalize
I1 = toGray01(A.I); I2 = toGray01(B.I);

% resize to common
H = min(size(I1,1), size(I2,1));
W = min(size(I1,2), size(I2,2));
I1 = imresize(I1,[H W]); I2 = imresize(I2,[H W]);

x1 = linspace(A.x(1), A.x(end), W);
x2 = linspace(B.x(1), B.x(end), W);
x  = (x1 + x2)/2;  % common x-axis in plot units

% normalize to [0,1] for screenshot-like stability
I1n = normalize01(I1);
I2n = normalize01(I2);

D = I1n - I2n;

out.mse = mean(D(:).^2);
out.psnr_db = 10*log10(1/(out.mse + 1e-12));
out.corr = corr(I1n(:), I2n(:));
out.mad = mean(abs(D(:)));
out.absdiff_p95 = prctile(abs(D(:)),95);
try
    out.ssim = ssim(I1n, I2n);
catch
    out.ssim = NaN;
end

% PSF center by brightest
[cy1,cx1] = ind2sub(size(I1n), argmax2(I1n));
[cy2,cx2] = ind2sub(size(I2n), argmax2(I2n));
cy = round((cy1+cy2)/2); cx = round((cx1+cx2)/2);

% FWHM lateral in axis units (same x units as image)
out.FWHM_lat_A = fwhm_1d_axis(I1n(cy,:), x, cx);
out.FWHM_lat_B = fwhm_1d_axis(I2n(cy,:), x, cx);

% PSLR/ISLR from -6 dB threshold region (>= 0.5 peak) connected to center
[out.PSLR_dB_A, out.ISLR_dB_A, mainMaskA] = pslr_islr_2d(I1n, cx, cy, 0.5);
[out.PSLR_dB_B, out.ISLR_dB_B, mainMaskB] = pslr_islr_2d(I2n, cx, cy, 0.5);

if showFigures
    out.fig_roi = figure('Name',sprintf('2D ROI %d',figIndex));
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile; imagesc(I1n); axis image off; colormap gray; title("A (yours)");
    hold on; plot(cx,cy,'r+'); hold off;
    nexttile; imagesc(I2n); axis image off; colormap gray; title("B (senior)");
    hold on; plot(cx,cy,'r+'); hold off;

    out.fig_overlay = figure('Name',sprintf('2D Overlay %d',figIndex));
    imagesc(alpha*I1n + (1-alpha)*I2n); axis image off; colormap gray;
    title(sprintf('Overlay alpha=%.2f',alpha));

    out.fig_absdiff = figure('Name',sprintf('2D AbsDiff %d',figIndex));
    imagesc(abs(D)); axis image off; colormap gray; title('|Me-MingChi|');

    out.fig_profiles = figure('Name',sprintf('2D Lateral profile %d',figIndex));
    plot(x, I1n(cy,:),'DisplayName','Me'); hold on;
    plot(x, I2n(cy,:),'DisplayName','MingChi'); grid on; legend('Location','best');
    title('Lateral profile at mainlobe depth (normalized)');

    out.fig_psf = figure('Name',sprintf('2D PSLR/ISLR masks %d',figIndex));
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile; imagesc(I1n); axis image off; colormap gray;
    hold on; contour(mainMaskA,[0.5 0.5],'r','LineWidth',1); plot(cx,cy,'r+'); hold off;
    title(sprintf('Me: FWHM=%.3g, PSLR=%.1fdB, ISLR=%.1fdB', out.FWHM_lat_A, out.PSLR_dB_A, out.ISLR_dB_A));
    nexttile; imagesc(I2n); axis image off; colormap gray;
    hold on; contour(mainMaskB,[0.5 0.5],'r','LineWidth',1); plot(cx,cy,'r+'); hold off;
    title(sprintf('MingChi: FWHM=%.3g, PSLR=%.1fdB, ISLR=%.1fdB', out.FWHM_lat_B, out.PSLR_dB_B, out.ISLR_dB_B));
else
    out.fig_roi=[]; out.fig_overlay=[]; out.fig_absdiff=[]; out.fig_profiles=[]; out.fig_psf=[];
end

out.kind = "image2d";
end

% ======================================================================
%                          1D LINE COMPARE
% ======================================================================

function out = compare_line1d(A, B, figIndex, showFigures)
% Beam pattern from your scripts are already in dB and peak-normalized (0 dB at max).
% We'll compute:
% - similarity on y (interpolated to common x)
% - FWHM in x units (width above -6 dB)
% - PSLR: highest sidelobe peak (in dB) outside mainlobe
% - ISLR: energy ratio sidelobe/mainlobe in linear domain

x1=A.x; y1=A.y;
x2=B.x; y2=B.y;

% common x range overlap
xmin = max(min(x1), min(x2));
xmax = min(max(x1), max(x2));
N = min(numel(x1), numel(x2));
x = linspace(xmin, xmax, max(200, N));

y1i = interp1(x1, y1, x, 'linear', 'extrap');
y2i = interp1(x2, y2, x, 'linear', 'extrap');

D = y1i - y2i;
out.mse = mean(D.^2);
out.psnr_db = 10*log10(1/(out.mse + 1e-12)); % not super meaningful for dB curves, but ok for comparison
out.corr = corr(y1i(:), y2i(:));
out.mad = mean(abs(D));
out.absdiff_p95 = prctile(abs(D),95);

% FWHM: width above -6 dB around peak (peak should be ~0 dB)
out.FWHM_A = fwhm_db_curve(x, y1i, -6);
out.FWHM_B = fwhm_db_curve(x, y2i, -6);

% PSLR / ISLR from 1D
[out.PSLR_dB_A, out.ISLR_dB_A, mainMaskA] = pslr_islr_1d_db(x, y1i, -6);
[out.PSLR_dB_B, out.ISLR_dB_B, mainMaskB] = pslr_islr_1d_db(x, y2i, -6);

if showFigures
    out.fig_roi = figure('Name',sprintf('1D curve %d',figIndex));
    plot(x,y1i,'DisplayName','Me'); hold on; plot(x,y2i,'DisplayName','MingChi'); grid on; legend;
    title(sprintf('1D beam pattern compare (Fig %d)',figIndex));

    out.fig_overlay = figure('Name',sprintf('1D overlay %d',figIndex));
    plot(x,y1i,'DisplayName','Me'); hold on;
    plot(x,y2i,'DisplayName','MingChi');
    yline(-6,'--','-6 dB');
    grid on; legend; title('Overlay with -6 dB line');

    out.fig_absdiff = figure('Name',sprintf('1D absdiff %d',figIndex));
    plot(x,abs(D)); grid on; title('|A-B| (dB)');

    out.fig_profiles = figure('Name',sprintf('1D mainlobe mask %d',figIndex));
    plot(x,y1i,'DisplayName','Me'); hold on;
    plot(x,y2i,'DisplayName','MingChi');
    scatter(x(mainMaskA), y1i(mainMaskA), 6, 'filled','DisplayName','Me mainlobe');
    scatter(x(mainMaskB), y2i(mainMaskB), 6, 'filled','DisplayName','MingChi mainlobe');
    grid on; legend('Location','best');
    title(sprintf('FWHM_A=%.3g, PSLR_A=%.1fdB, ISLR_A=%.1fdB', out.FWHM_A, out.PSLR_dB_A, out.ISLR_dB_A));

    out.fig_psf = figure('Name',sprintf('1D summary %d',figIndex));
    plot(x,y1i,'DisplayName','Me'); hold on; plot(x,y2i,'DisplayName','MingChi'); grid on; legend;
    title(sprintf('A: FWHM=%.3g PSLR=%.1f ISLR=%.1f | B: FWHM=%.3g PSLR=%.1f ISLR=%.1f', ...
        out.FWHM_A,out.PSLR_dB_A,out.ISLR_dB_A, out.FWHM_B,out.PSLR_dB_B,out.ISLR_dB_B));
else
    out.fig_roi=[]; out.fig_overlay=[]; out.fig_absdiff=[]; out.fig_profiles=[]; out.fig_psf=[];
end

out.kind = "line1d";
end

% ======================================================================
%                               SAVE
% ======================================================================

function save_all_figs(out, saveDir, baseName)
if ~isempty(out.fig_roi),     saveas(out.fig_roi,     fullfile(saveDir, baseName+"_roi.png")); end
if ~isempty(out.fig_overlay), saveas(out.fig_overlay, fullfile(saveDir, baseName+"_overlay.png")); end
if ~isempty(out.fig_absdiff), saveas(out.fig_absdiff, fullfile(saveDir, baseName+"_absdiff.png")); end
if ~isempty(out.fig_profiles),saveas(out.fig_profiles,fullfile(saveDir, baseName+"_profiles.png")); end
if ~isempty(out.fig_psf),     saveas(out.fig_psf,     fullfile(saveDir, baseName+"_psf.png")); end
end

% ======================================================================
%                          METRIC HELPERS
% ======================================================================

function I = toGray01(I)
I = im2double(I);
if ndims(I)==3, I = rgb2gray(I); end
end

function R = normalize01(A)
A = double(A);
R = (A - min(A(:))) / (max(A(:)) - min(A(:)) + eps);
end

function idx = argmax2(M)
[~,idx] = max(M(:));
end

function w = fwhm_1d_axis(profile, x, centerIdx)
p = double(profile(:))';
pk = p(centerIdx);
if pk <= 0, w = NaN; return; end
th = 0.5*pk;

l = centerIdx;
while l>1 && p(l) >= th, l = l-1; end
r = centerIdx;
while r<numel(p) && p(r) >= th, r = r+1; end

w = x(r) - x(l);
end

function [PSLR_dB, ISLR_dB, mainMask] = pslr_islr_2d(I, cx, cy, frac)
pk = I(cy,cx);
th = frac*pk;

BW = I >= th;

% main connected component that contains center pixel
CC = bwconncomp(BW, 8);
mainMask = false(size(BW));
if CC.NumObjects == 0
    PSLR_dB = NaN; ISLR_dB = NaN; return;
end

centerIdx = sub2ind(size(BW), cy, cx);
found = false;
for k=1:CC.NumObjects
    if any(CC.PixelIdxList{k} == centerIdx)
        mainMask(CC.PixelIdxList{k}) = true;
        found = true;
        break;
    end
end
if ~found
    mainMask = BW;
end

mainPeak = max(I(mainMask), [], 'all');
sidePeak = max(I(~mainMask), [], 'all');

PSLR_dB = 20*log10((sidePeak + eps) / (mainPeak + eps));

Emain = sum(I(mainMask).^2, 'all');
Eside = sum(I(~mainMask).^2, 'all');
ISLR_dB = 10*log10((Eside + eps) / (Emain + eps));
end

function w = fwhm_db_curve(x, y_db, level_db)
% width of region >= level_db around global max
[~,i0] = max(y_db);
l=i0; while l>1 && y_db(l) >= level_db, l=l-1; end
r=i0; while r<numel(y_db) && y_db(r) >= level_db, r=r+1; end
w = x(r) - x(l);
end

function [PSLR_dB, ISLR_dB, mainMask] = pslr_islr_1d_db(x, y_db, level_db)
% mainMask: contiguous region around peak where y >= level_db
[~,i0] = max(y_db);
mainMask = false(size(y_db));
l=i0; while l>1 && y_db(l) >= level_db, l=l-1; end
r=i0; while r<numel(y_db) && y_db(r) >= level_db, r=r+1; end
mainMask(l:r) = true;

% PSLR: max sidelobe peak (in linear amplitude ratio, expressed in dB)
% since y_db is already relative to main peak (~0 dB), PSLR is simply max(y_db outside mainlobe)
PSLR_dB = max(y_db(~mainMask));

% ISLR: energy ratio (linear)
A = 10.^(y_db/20); % amplitude
Emain = sum(A(mainMask).^2);
Eside = sum(A(~mainMask).^2);
ISLR_dB = 10*log10((Eside + eps)/(Emain + eps));
end

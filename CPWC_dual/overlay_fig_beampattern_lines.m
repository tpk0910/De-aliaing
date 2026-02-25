function overlay_fig_beampattern_lines(figFiles, labels, doNormalize)
% overlay_fig_beampattern_lines
% ------------------------------------------------------------
% Overlay 1D beam pattern curves saved in .fig files into one plot.
%
% figFiles    : cell array of .fig file paths
% labels      : cell array of legend labels (same length as figFiles)
% doNormalize : true -> normalize each curve to its own max (0 dB)
%
% Usage:
% figFiles = {'case48.fig','case64.fig','case96.fig'};
% labels   = {'48 el','64 el','96 el'};
% overlay_fig_beampattern_lines(figFiles, labels, true);

if nargin < 2 || isempty(labels)
    labels = figFiles;
end
if nargin < 3
    doNormalize = true;
end

assert(numel(figFiles) == numel(labels), 'labels must match figFiles length.');

% New figure for overlay
fOut = figure('Name','Beam patterns','Color','w');
axOut = axes('Parent', fOut); hold(axOut, 'on'); grid(axOut, 'on');

allX = [];
allY = [];

for i = 1:numel(figFiles)
    % Open .fig invisibly
    f = openfig(figFiles{i}, 'invisible');

    % Find line objects (beam patterns usually are lines)
    ax = findobj(f, 'Type', 'axes');
    ln = findobj(ax, 'Type', 'line');

    if isempty(ln)
        close(f);
        error('No line objects found in %s. If this is an image figure, use the 2D method.', figFiles{i});
    end

    % Many figs may have multiple lines; usually the "main" one has most points
    [~, idx] = max(arrayfun(@(h) numel(h.XData), ln));
    mainLine = ln(idx);

    x = mainLine.XData(:);
    y = mainLine.YData(:);

    % Optional normalize to 0 dB peak for fair comparison
    if doNormalize
        y = y - max(y);  % assumes y is already in dB
    end

    plot(axOut, x, y, 'LineWidth', 1.5, 'DisplayName', labels{i});

    allX = [allX; x];
    allY = [allY; y];

    close(f);
end

xlabel(axOut, 'Angle / Lateral');
ylabel(axOut, 'Amplitude (dB)');
legend(axOut, 'Location', 'best');
title(axOut, 'Beam pattern');

% Set reasonable axes limits if data exist
if ~isempty(allX)
    xlim(axOut, [min(allX) max(allX)]);
end
if ~isempty(allY)
    ylim(axOut, [max(min(allY), -80) 5]); % clip to [-80, 5] dB range
end

end
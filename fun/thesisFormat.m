function thesisFormat(fontName, subH, fig)
%THESISFORMAT Apply consistent thesis-ready typography and figure sizing.
%   THESISFORMAT formats the current figure. It preserves axes limits,
%   direction, label/title text, and data-specific line colors and styles.
%   THESISFORMAT(FONTNAME, SUBH) sets the font and approximate subplot height
%   in centimeters. Figures with overlaid axes are counted as one position.
%   THESISFORMAT(FONTNAME, SUBH, FIG) formats an explicit figure handle.
%   Pass [] for either optional style argument to use its default.
%   Figure dimensions (including SUBH) use a 1.25 scale.
%   Default single-panel height after scaling is 6 cm.
%   Main axes and labels use 14 pt; inset axes use 10 pt.
%   Colorbars use 16 pt; legends use 12 pt.
%   Insets and their tick labels are kept inside the main axes plot box.
%   To associate an inset explicitly (including one initially outside), use
%   setappdata(insetAxes, 'thesisFormatMainAxes', mainAxes). Both axes must
%   share a parent. Otherwise, smaller contained axes are detected as insets.
%   Base subplot width/height ratio is 1.75.
%   Axes boxes are on; legends have white boxes with light gray borders.
%   Uses inward major/minor ticks and light gray solid major grids.
%   Minor grid lines are disabled.
%   Call after creating all axes, labels, titles, and legends.

formatScale = 1.25;

if nargin < 1 || isempty(fontName)
    fontName = 'Times New Roman';
end
if nargin < 2 || isempty(subH)
    subH = 4.8;
end

if nargin < 3 || isempty(fig)
    fig = gcf;
end
allAxes = findall(fig, 'Type', 'axes');
if isempty(allAxes)
    return
end

% Exclude axes-like containers such as legends and colorbars.
keep = true(size(allAxes));
for k = 1:numel(allAxes)
    keep(k) = ~strcmp(allAxes(k).Tag, 'legend') ...
        && ~strcmp(allAxes(k).Tag, 'Colorbar');
end
allAxes = allAxes(keep);
if isempty(allAxes)
    return
end

% Compare normalized positions so physical axes units do not affect sizing.
for k = 1:numel(allAxes)
    allAxes(k).Units = 'normalized';
end
% Insets are styled but do not count as additional subplot rows/columns.
% Coincident positions represent overlaid axes, not extra panels.
positions = vertcat(allAxes.Position);
isInset = false(size(allAxes));
insetMainIndex = zeros(size(allAxes));
for k = 1:numel(allAxes)
    mainAxes = getappdata(allAxes(k), 'thesisFormatMainAxes');
    if isscalar(mainAxes) && isgraphics(mainAxes, 'axes')
        mainIndex = find(allAxes == mainAxes, 1);
        if ~isempty(mainIndex) && mainIndex ~= k ...
                && isequal(allAxes(k).Parent, mainAxes.Parent)
            isInset(k) = true;
            insetMainIndex(k) = mainIndex;
            continue
        end
    end
    pos = positions(k, :);
    contained = positions(:, 1) <= pos(1) + 1e-3 ...
        & positions(:, 2) <= pos(2) + 1e-3 ...
        & positions(:, 1) + positions(:, 3) >= pos(1) + pos(3) - 1e-3 ...
        & positions(:, 2) + positions(:, 4) >= pos(2) + pos(4) - 1e-3;
    larger = positions(:, 3) > pos(3) + 1e-3 ...
        & positions(:, 4) > pos(4) + 1e-3;
    sameParent = arrayfun(@(ax) isequal(ax.Parent, allAxes(k).Parent), allAxes);
    candidates = find(contained & larger & sameParent);
    if ~isempty(candidates)
        % Associate with the smallest enclosing axes for nested layouts.
        [~, nearest] = min(prod(positions(candidates, 3:4), 2));
        isInset(k) = true;
        insetMainIndex(k) = candidates(nearest);
    end
end
% Merge small layout differences when counting distinct subplot columns/rows.
lefts = unique(round(positions(~isInset, 1), 3));
bases = unique(round(positions(~isInset, 2), 3));

subplotAspectRatio = 1.75; % Base subplot width/height ratio.
subW = subplotAspectRatio * subH;
totalW = subW * numel(lefts);
totalH = subH * numel(bases);

% Overlaid top and bottom x-axes need extra physical height for both sets
% of tick labels and axis labels, especially when the figure is exported.
hasTopXAxis = false;
for k = 1:numel(allAxes)
    if strcmp(allAxes(k).XAxisLocation, 'top')
        hasTopXAxis = true;
        break
    end
end
if hasTopXAxis
    totalH = max(totalH, 8);
end
totalW = formatScale * totalW;
totalH = formatScale * totalH;

fig.Units = 'centimeters';
fig.Position = [10 3 totalW totalH];
fig.Color = [1 1 1];
fig.PaperUnits = 'centimeters';
fig.PaperPositionMode = 'auto';
fig.PaperSize = [totalW totalH];

% Reserve room above and below axes that have a top x-axis. A larger figure
% alone does not prevent labels from clipping when axes fill the figure.
originalPositions = vertcat(allAxes.Position);
for k = 1:numel(allAxes)
    if ~strcmp(allAxes(k).XAxisLocation, 'top')
        continue
    end
    pos = originalPositions(k, :);
    % Move coincident axes together so overlaid coordinate systems stay aligned.
    sameAxes = all(abs(originalPositions - pos) < 1e-3, 2);
    lowerEdge = max(pos(2), 0.16);
    upperEdge = min(pos(2) + pos(4), 0.74);
    if upperEdge - lowerEdge >= 0.4
        adjustedPos = [pos(1), lowerEdge, pos(3), upperEdge-lowerEdge];
        for j = find(sameAxes).'
            allAxes(j).Position = adjustedPos;
        end
    end
end

rulerNames = {'XAxis', 'YAxis', 'ZAxis'};
for k = 1:numel(allAxes)
    ax = allAxes(k);
    axisFontSize = 14;
    if isInset(k)
        axisFontSize = 10;
    end
    ax.FontName = fontName;
    ax.FontUnits = 'points';
    ax.FontSize = axisFontSize;
    ax.FontWeight = 'normal';
    ax.Box = 'on';
    ax.LineWidth = 0.7;
    ax.TickLabelInterpreter = 'latex';
    ax.Color = 'white';
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.ZColor = 'black';
    ax.TickDir = 'in';
    ax.TickLength = [0.015, 0.025];
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.ZMinorTick = 'on';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.ZGrid = 'on';
    ax.GridLineStyle = '-';
    ax.GridColor = [0.85, 0.85, 0.85];
    ax.GridAlpha = 1;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    ax.ZMinorGrid = 'off';
    ax.TitleFontWeight = 'normal';
    ax.Title.FontName = fontName;
    ax.Title.Interpreter = 'latex';
    if isInset(k)
        ax.Title.FontSize = axisFontSize;
    end

    % Style labels without changing their text or the axis configuration.
    for r = 1:numel(rulerNames)
        if isprop(ax, rulerNames{r})
            rulers = ax.(rulerNames{r});
            for j = 1:numel(rulers)
                if isprop(rulers(j), 'FontSize')
                    rulers(j).FontSize = axisFontSize;
                end
                if isprop(rulers(j), 'Label')
                    rulers(j).Label.FontName = fontName;
                    rulers(j).Label.FontUnits = 'points';
                    rulers(j).Label.FontSize = axisFontSize;
                    rulers(j).Label.Interpreter = 'latex';
                end
            end
        end
    end

    % Increase default-width line plots while preserving explicitly chosen
    % widths and leaving limits, title text, and colors untouched.
    lines = findall(ax, 'Type', 'line');
    for j = 1:numel(lines)
        if lines(j).LineWidth == 1
            lines(j).LineWidth = 2;
        end
    end
end

% Format legends and colorbars when present without changing their location.
legends = findall(fig, 'Type', 'legend');
for k = 1:numel(legends)
    legends(k).FontName = fontName;
    legends(k).FontUnits = 'points';
    legends(k).FontSize = 12;
    legends(k).Interpreter = 'latex';
    legends(k).ItemTokenSize = [15, 18]; % Line sample dimensions in points.
    legends(k).Box = 'on';
    legends(k).Color = 'white';
    legends(k).EdgeColor = [0.8, 0.8, 0.8];
    legends(k).LineWidth = 0.5;
end

colorbars = findall(fig, 'Type', 'colorbar');
for k = 1:numel(colorbars)
    colorbars(k).FontName = fontName;
    colorbars(k).FontUnits = 'points';
    colorbars(k).FontSize = 16;
    colorbars(k).Label.FontName = fontName;
    colorbars(k).Label.FontUnits = 'points';
    colorbars(k).Label.FontSize = 16;
    colorbars(k).Label.Interpreter = 'latex';
end

% Render pending layout changes before reading text extents through TightInset.
drawnow;
if ~any(isInset)
    return
end
% TightInset includes tick labels and axis labels. Measure after rendering
% at the final figure size, then constrain their full extent to the main box.
% A second pass accounts for tick-label changes caused by resizing an inset.
for pass = 1:2
    for k = find(isInset).'
        ax = allAxes(k);
        mainPos = allAxes(insetMainIndex(k)).Position;
        padding = 0.02 * mainPos(3:4);
        tight = ax.TightInset;
        lower = mainPos(1:2) + padding + tight(1:2);
        upper = mainPos(1:2) + mainPos(3:4) - padding - tight(3:4);
        pos = ax.Position;
        pos(3:4) = min(pos(3:4), max(upper - lower, eps));
        pos(1:2) = max(lower, min(pos(1:2), upper - pos(3:4)));
        ax.PositionConstraint = 'innerposition';
        ax.Position = pos;
    end
    drawnow;
end
end

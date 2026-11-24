function stateFigure = plotStateOptimization(figureName, cellSOC, stateOCV, ...
        absoluteEndpoints, relativeEndpoints, ocpData, gemColors)
%PLOTSTATEOPTIMIZATION Compare lithiation coordinates for one fitted state.
%   fig = plotStateOptimization(name, cellSOC, stateOCV, absoluteEndpoints, ...
%       relativeEndpoints, ocpData, gemColors) creates a three-panel figure.
%   cellSOC and stateOCV [V] are matching vectors. Endpoint arrays are 2-by-2:
%   rows [SOC0; SOC1], columns [negative positive], dimensionless lithiation.
%   ocpData has positiveThetaAbsolute, positiveThetaRelative, positiveOCP,
%   negativeThetaAbsolute, negativeThetaRelative, negativeOCP, lithiationLimits
%   and socLimits. OCP vectors are in V; coordinates and two-element limits
%   are dimensionless. Each electrode's three LUT vectors have matching sizes.
%   gemColors is an N-by-3 RGB palette with N >= 5 and values in [0,1].
%   Returns the formatted three-panel figure handle.
stateFigure = figure('Name', figureName);
stateAxes = gobjects(3, 1);
for iPanel = 1:3
    stateAxes(iPanel) = subplot(3, 1, iPanel, 'Parent', stateFigure);
    hold(stateAxes(iPanel), 'on');
end
plot(stateAxes(1), ocpData.positiveThetaAbsolute, ocpData.positiveOCP, '-', ...
    'Color', gemColors(4,:), 'LineWidth', 2, 'DisplayName', 'Abs. Lit.');
plot(stateAxes(1), ocpData.positiveThetaRelative, ocpData.positiveOCP, '--', ...
    'Color', gemColors(4,:), 'LineWidth', 2, 'DisplayName', 'Rel. Lit.');
ylabel(stateAxes(1), '$U_\mathrm{ocp}^\mathrm{p}$ (V)', 'Interpreter', 'latex');
legend(stateAxes(1), 'Location', 'best', 'Interpreter', 'latex');
plot(stateAxes(2), cellSOC, stateOCV, ...
    'Color', gemColors(1,:), 'LineWidth', 2);
ylabel(stateAxes(2), '$U_\mathrm{ocv}$ (V)', 'Interpreter', 'latex');
plot(stateAxes(3), ocpData.negativeThetaAbsolute, ocpData.negativeOCP, '-', ...
    'Color', gemColors(5,:), 'LineWidth', 2, 'DisplayName', 'Abs. Lit.');
plot(stateAxes(3), ocpData.negativeThetaRelative, ocpData.negativeOCP, '--', ...
    'Color', gemColors(5,:), 'LineWidth', 2, 'DisplayName', 'Rel. Lit.');
ylabel(stateAxes(3), '$U_\mathrm{ocp}^\mathrm{n}$ (V)', 'Interpreter', 'latex');
% Match panel order [positive OCP; cell OCV; negative OCP], keeping columns
% in SOC endpoint order even though positive stoichiometry decreases on charge.
absoluteGuideLocations = [absoluteEndpoints(:,2).'; 0, 1; absoluteEndpoints(:,1).'];
relativeGuideLocations = [relativeEndpoints(:,2).'; ...
    0, 1; relativeEndpoints(:,1).'];
for iPanel = 1:3
    if iPanel == 2
        xlim(stateAxes(iPanel), ocpData.socLimits);
        xlabel(stateAxes(iPanel), 'SOC z');
    else
        xlim(stateAxes(iPanel), ocpData.lithiationLimits);
        xlabel(stateAxes(iPanel), 'Lithiation $\theta$, $\tilde{\theta}$', ...
            'Interpreter', 'latex');
    end
    guideLocations = absoluteGuideLocations(iPanel,:);
    guideStyles = {':'};
    guideColors = [0.65, 0.65, 0.65];
    if iPanel ~= 2
        guideLocations = [guideLocations; relativeGuideLocations(iPanel,:)];
        guideStyles = {'-', '--'};
        guideColors = [0, 0, 0; 0.65, 0.65, 0.65];
    end
    for iCoordinate = 1:size(guideLocations, 1)
        for iEndpoint = 1:2
            % Extend labels toward the nearest edge to keep the central area clear.
            guideHorizontalAlignment = 'right';
            if guideLocations(iCoordinate,iEndpoint) > mean(xlim(stateAxes(iPanel)))
                guideHorizontalAlignment = 'left';
            end
            guideVerticalAlignment = 'bottom';
            if iPanel == 1
                guideVerticalAlignment = 'top';
            end
            xline(stateAxes(iPanel), guideLocations(iCoordinate,iEndpoint), ...
                guideStyles{iCoordinate}, sprintf('$z=%d$', iEndpoint-1), ...
                'Color', guideColors(iCoordinate,:), ...
                'LineWidth', 1, 'HandleVisibility', 'off', 'Interpreter', 'latex', ...
                'LabelOrientation', 'horizontal', ...
                'LabelHorizontalAlignment', guideHorizontalAlignment, ...
                'LabelVerticalAlignment', guideVerticalAlignment);
        end
    end
end
thesisFormat([], [], stateFigure);
end

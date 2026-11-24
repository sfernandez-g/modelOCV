%MODELOCV Compare absolute and relative LLI for synthetic OCV states.
%   Run modelOCV.m to display tables, open five figures, and retain results.
%   Requires MATLAB R2023b or newer and Optimization Toolbox.
%   Clears the workspace and closes existing figures.
%   Rows/columns of thetaLimits(:,:,iState): [SOC0; SOC1] / [negative positive].
%   reconstructedOCV(:,iState) contains OCV [V] on the shared cellSOC grid.
%   Generates the Section 3.5 table and capacity/SOC plots, including insets.
%   results retains full-precision states, calibration, and solver diagnostics.
%   Displays tables in the Command Window and generates figures without saving files.

clc; clear; close all

% Resolve project paths from this file so execution does not depend on pwd.
scriptDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDirectory, 'fun'));
dataDirectory = fullfile(scriptDirectory, 'input');
gemColors = orderedcolors('gem');

%% Load electrode data and define OCP functions
% Use only the supplied numerical average-branch LUTs, not stored handles.
negativeReference = load(fullfile(dataDirectory, 'Uocp_n.mat'), 'model');
positiveReference = load(fullfile(dataDirectory, 'Uocp_p.mat'), 'model');
negativeElectrodeOCP = negativeReference.model.LUT.U_avg(:);
negativeThetaAbs = negativeReference.model.LUT.theta_avg(:);
negativeThetaRel = negativeReference.model.LUT.thetarel_avg(:);
positiveElectrodeOCP = positiveReference.model.LUT.U_avg(:);
positiveThetaAbs = positiveReference.model.LUT.theta_avg(:);
positiveThetaRel = positiveReference.model.LUT.thetarel_avg(:);
referenceProvenance = struct('negativeSource', 'Uocp_n.mat', ...
    'positiveSource', 'Uocp_p.mat', ...
    'arrays', 'model.LUT.theta_avg, thetarel_avg, U_avg', ...
    'conversion', 'Unmodified average-branch arrays; no fitting or resampling.');
clear negativeReference positiveReference

% Fixed half-cell calibration: thetaAbs = a + b*thetaRel for all states.
% Use the supplied LUT pairs rather than normalizing each operating window.
[negativeOrigin, negativeScale] = coordinateCalibration(negativeThetaAbs, negativeThetaRel);
[positiveOrigin, positiveScale] = coordinateCalibration(positiveThetaAbs, positiveThetaRel);
coordinateOrigin = [negativeOrigin, positiveOrigin];
coordinateScale = [negativeScale, positiveScale];
coordinateCalibrationSummary = table(["Negative"; "Positive"], ...
    coordinateOrigin.', coordinateScale.', ...
    'VariableNames', {'Electrode', 'AbsoluteOrigin', 'AbsoluteScale'});

negativeOCPfun = prepareOCPInterpolant(negativeThetaAbs, negativeElectrodeOCP);
positiveOCPfun = prepareOCPInterpolant(positiveThetaAbs, positiveElectrodeOCP);

%% Fresh reference and illustrative degradation scenarios
cellSOC = linspace(0, 1, 1e3);
cutOffVoltage = [2.5, 4.2];       % [minimum maximum], V
negativeElectrodeCapacity = 3.5;    % fresh full-stoichiometry capacity, Ah
positiveElectrodeCapacity = 4;  % fresh full-stoichiometry capacity, Ah
thetaGuess = [0.1, 0.9; 0.9, 0.1]; % fresh inventory reference and initial guess
% Resolve cutoff residuals below 1 microvolt while retaining exact balances.
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', ...
    'FiniteDifferenceType', 'central', 'FiniteDifferenceStepSize', 1e-8, ...
    'OptimalityTolerance', 1e-9, 'ConstraintTolerance', 1e-10, ...
    'StepTolerance', 1e-12, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 10000);

% Columns [LLI LAMn LAMp]; each aged row is relative to the fresh state.
% These are assumed comparison scenarios, not measured degradation values.
% Draw two mixed-loss scenarios uniformly from the 5--20% grid.
randomSeed = 2027;
randomLossRange = [0.05, 0.20]; % Fractions of the fresh reference.
randomLossStep = 0.05;
randomLossLevels = randomLossRange(1):randomLossStep:randomLossRange(2);
% A local stream makes the scenarios repeatable without resetting MATLAB's RNG.
scenarioStream = RandStream('mt19937ar', 'Seed', randomSeed);
mixedLosses = randomLossLevels(randi(scenarioStream, numel(randomLossLevels), 2, 3));
degradationModes = [0,    0,    0;
                    0.05, 0,    0;
                    0,    0,    0.10;
                    mixedLosses];
nStates = size(degradationModes, 1);
stateNames = ["f"; compose("a%d", (1:nStates-1).')];
thetaLimits = zeros(2, 2, nStates);
thetaGuesses = zeros(2, 2, nStates);
reconstructedOCV = zeros(numel(cellSOC), nStates);
cellCapacity = zeros(nStates, 1);
electrodeCapacity = zeros(nStates, 2);
exitflag = zeros(nStates, 1);
residuals = cell(nStates, 1);
maxCutoffError = zeros(nStates, 1); % V
capacityMismatch = zeros(nStates, 1); % Ah
inventoryMismatch = zeros(nStates, 1); % Ah
extrapolationSpan = zeros(nStates, 2); % Maximum excursion outside each LUT [-].
voltageTolerance = 1e-6; % V; diagnostic threshold, does not alter the solver.
capacityTolerance = 1e-6*min(negativeElectrodeCapacity, positiveElectrodeCapacity); % Ah

%% Build the fresh state first, then each aged state independently
for iState = 1:nStates
    if iState == 1
        stateGuess = thetaGuess;
    else
        stateGuess = freshThetaLimits;
    end
    thetaGuesses(:,:,iState) = stateGuess;
    [stateTheta, stateOCV, stateCapacity, stateElectrodeCapacity, ...
        stateExitflag, stateResidual] = buildOCV(negativeOCPfun, positiveOCPfun, ...
        cutOffVoltage, negativeElectrodeCapacity, positiveElectrodeCapacity, ...
        stateGuess, cellSOC, degradationModes(iState,1), ...
        degradationModes(iState,2), degradationModes(iState,3), options);
    if stateExitflag <= 0 || ~isfinite(stateCapacity) || stateCapacity <= 0 ...
            || any(~isfinite(stateTheta(:))) || any(~isfinite(stateOCV(:))) ...
            || stateTheta(2,1) <= stateTheta(1,1) ...
            || stateTheta(1,2) <= stateTheta(2,2)
        error('modelOCV:BuildFailed', '%s failed (exitflag %d): %s', ...
            stateNames(iState), stateExitflag, stateResidual.optimization.message);
    end
    thetaLimits(:,:,iState) = stateTheta;
    reconstructedOCV(:,iState) = stateOCV(:);
    cellCapacity(iState) = stateCapacity;
    electrodeCapacity(iState,:) = stateElectrodeCapacity;
    exitflag(iState) = stateExitflag;
    residuals{iState} = stateResidual;
    maxCutoffError(iState) = max(abs([stateResidual.voltageEoD, stateResidual.voltageEoC]));
    capacityMismatch(iState) = stateResidual.capacityMismatch;
    inventoryMismatch(iState) = stateResidual.lithiumInventory;
    measuredLimits = [min(negativeThetaAbs), min(positiveThetaAbs); ...
        max(negativeThetaAbs), max(positiveThetaAbs)];
    % Keep the largest endpoint excursion for each electrode; zero means none.
    extrapolationSpan(iState,:) = max([zeros(1, 2); ...
        measuredLimits(1,:)-min(stateTheta, [], 1); ...
        max(stateTheta, [], 1)-measuredLimits(2,:)], [], 1);
    if maxCutoffError(iState) > voltageTolerance ...
            || abs(capacityMismatch(iState)) > capacityTolerance ...
            || abs(inventoryMismatch(iState)) > capacityTolerance
        error('modelOCV:ResidualTolerance', ...
            '%s: cutoff error %.3g V, charge mismatch %.3g Ah, inventory mismatch %.3g Ah.', ...
            stateNames(iState), maxCutoffError(iState), ...
            capacityMismatch(iState), inventoryMismatch(iState));
    end
    if iState == 1
        % This fitted absolute reference remains fixed for every aged call.
        freshThetaLimits = stateTheta;
    end
end

%% Coordinate analysis: same fitted OCV states, without refitting
% theta arrays: [SOC endpoint, electrode, state]; capacities: [state, electrode].
% All inventories below are charge equivalents in Ah, not amounts in mol.
% Shape the calibration as [1, electrode, 1] to apply it across endpoints/states.
thetaLimitsRelative = (thetaLimits-reshape(coordinateOrigin, 1, 2, 1)) ...
    ./reshape(coordinateScale, 1, 2, 1);
% Fixed affine scaling gives Qtilde = b*Qabsolute (PDF Eq. 22).
relativeElectrodeCapacity = electrodeCapacity.*coordinateScale;
inventoryAbsoluteEndpoints = zeros(nStates, 2);
inventoryRelativeEndpoints = zeros(nStates, 2);
for iEndpoint = 1:2
    % Turn each endpoint slice into [state, electrode] to match the capacities.
    absoluteEndpoints = reshape(thetaLimits(iEndpoint,:,:), 2, nStates).';
    relativeEndpoints = reshape(thetaLimitsRelative(iEndpoint,:,:), 2, nStates).';
    inventoryAbsoluteEndpoints(:,iEndpoint) = sum(electrodeCapacity.*absoluteEndpoints, 2);
    inventoryRelativeEndpoints(:,iEndpoint) = sum(relativeElectrodeCapacity.*relativeEndpoints, 2);
end
% Origin inventory varies with LAM even though the calibration stays fixed.
inventoryOrigin = sum(electrodeCapacity.*coordinateOrigin, 2);
inventoryCorrectedEndpoints = inventoryRelativeEndpoints+inventoryOrigin;
% Evaluate at SOC=1 as in the paper; retain both endpoints for inspection.
% buildOCV constrains SOC=0; the small endpoint difference is solver residual.
inventoryAbsolute = inventoryAbsoluteEndpoints(:,2);
inventoryRelative = inventoryRelativeEndpoints(:,2);
inventoryCorrected = inventoryCorrectedEndpoints(:,2);
lliAbsolute = 1-inventoryAbsolute/inventoryAbsolute(1);
lliRelativeCorrected = 1-inventoryCorrected/inventoryCorrected(1);
% Scale the roundoff threshold by the uncanceled electrode contributions
% before dividing by a potentially small relative inventory.
relativeInventoryTolerance = 100*eps(max(1, sum(abs( ...
    relativeElectrodeCapacity(1,:).*reshape(thetaLimitsRelative(2,:,1), 1, 2)))));
if inventoryRelative(1) > relativeInventoryTolerance
    lliRelativeApparent = 1-inventoryRelative/inventoryRelative(1);
    % PDF Eq. 27: separate true inventory loss and capacity-dependent origins.
    lliRelativeFromBalance = ((inventoryAbsolute(1)-inventoryAbsolute) ...
        -sum((electrodeCapacity(1,:)-electrodeCapacity).*coordinateOrigin, 2)) ...
        /inventoryRelative(1);
else
    error('modelOCV:InvalidRelativeInventory', ...
        'Fresh relative inventory must be positive to report apparent fractional LLI.');
end
lamAbsolute = 1-electrodeCapacity./electrodeCapacity(1,:);
lamRelative = 1-relativeElectrodeCapacity./relativeElectrodeCapacity(1,:);
lliComparison = table(string(stateNames), 100*degradationModes(:,1), ...
    inventoryAbsolute, inventoryOrigin, inventoryRelative, inventoryCorrected, ...
    inventoryAbsolute(1)-inventoryAbsolute, inventoryRelative(1)-inventoryRelative, ...
    100*lliAbsolute, 100*lliRelativeApparent, 100*lliRelativeCorrected, ...
    100*(lliRelativeApparent-lliAbsolute), 100*lliRelativeFromBalance, ...
    'VariableNames', {'State', 'LLI_input_pct', 'InventoryAbsolute_Ah', ...
    'OriginInventory_Ah', 'InventoryRelative_Ah', 'InventoryCorrected_Ah', ...
    'LossAbsolute_Ah', 'LossRelative_Ah', 'LLI_absolute_pct', ...
    'LLI_relative_apparent_pct', 'LLI_relative_corrected_pct', ...
    'ApparentBias_pp', 'LLI_relative_balance_pct'});
lamComparison = table(string(stateNames), 100*lamAbsolute(:,1), ...
    100*lamRelative(:,1), 100*lamAbsolute(:,2), 100*lamRelative(:,2), ...
    'VariableNames', {'State', 'LAMn_absolute_pct', 'LAMn_relative_pct', ...
    'LAMp_absolute_pct', 'LAMp_relative_pct'});

%% Summarize capacities, degradation inputs, and optimization diagnostics
capacityRetention = 100*cellCapacity/cellCapacity(1);
% The outer product gives [SOC sample, state] discharged capacities.
dischargedCapacity = (1-cellSOC(:))*cellCapacity.'; % Ah since SOC=1
stateSummary = table(string(stateNames), 100*degradationModes(:,1), ...
    100*degradationModes(:,2), 100*degradationModes(:,3), cellCapacity, ...
    capacityRetention, electrodeCapacity(:,1), electrodeCapacity(:,2), ...
    1e3*maxCutoffError, capacityMismatch, inventoryMismatch, exitflag, ...
    'VariableNames', {'State', 'LLI_pct', 'LAMn_pct', 'LAMp_pct', ...
    'Capacity_Ah', 'Retention_pct', 'Qn_Ah', 'Qp_Ah', ...
    'MaxCutoffError_mV', 'ChargeMismatch_Ah', 'InventoryMismatch_Ah', 'Exitflag'});
%% One row per state, with the fields used in the manuscript table
negativeWindows = compose("%.3f--%.3f", reshape(thetaLimits(1,1,:), [], 1), ...
    reshape(thetaLimits(2,1,:), [], 1));
positiveWindows = compose("%.3f--%.3f", reshape(thetaLimits(1,2,:), [], 1), ...
    reshape(thetaLimits(2,2,:), [], 1));
paperTable = table(stateNames, cellCapacity, electrodeCapacity(:,1), ...
    electrodeCapacity(:,2), negativeWindows, positiveWindows, ...
    inventoryAbsolute, inventoryRelative, 100*lliAbsolute, ...
    100*lliRelativeApparent, 'VariableNames', {'State', 'Qc_Ah', 'Qn_Ah', ...
    'Qp_Ah', 'x0_x100', 'y0_y100', 'QLi_Ah', 'QLi_relative_Ah', ...
    'LLI_pct', 'LLI_relative_pct'});
disp(paperTable);
disp(stateSummary(:, {'State', 'LLI_pct', 'LAMn_pct', 'LAMp_pct', 'Retention_pct'}));
disp(coordinateCalibrationSummary);
extrapolationSummary = table(stateNames, extrapolationSpan(:,1), ...
    extrapolationSpan(:,2), 'VariableNames', ...
    {'State', 'NegativeThetaExcursion', 'PositiveThetaExcursion'});
disp(extrapolationSummary);

%% Electrode OCP data
figures = gobjects(5, 1);
figures(1) = figure('Name', 'Electrode OCP data');
ax = axes('Parent', figures(1));
hold(ax, 'on');
plot(ax, negativeThetaAbs, negativeElectrodeOCP, 'Color', gemColors(5,:), ...
    'LineWidth', 2, 'DisplayName', '$U_\mathrm{ocp}^\mathrm{n}$');
plot(ax, positiveThetaAbs, positiveElectrodeOCP, 'Color', gemColors(4,:), ...
    'LineWidth', 2, 'DisplayName', '$U_\mathrm{ocp}^\mathrm{p}$');
xlabel(ax, 'Absolute stoichiometry $\theta$');
ylabel(ax, 'Potential vs. Li (V)');
legend(ax, 'Location', 'best');
grid(ax, 'on');
thesisFormat([], [], figures(1));

%% Fresh optimization: compare absolute and relative lithiation, then cell SOC
% Electrode panels compare the numerical lithiation coordinates directly.
% Their guides mark SOC endpoints in both lithiation coordinates.
% The middle panel uses full-cell SOC; electrode curves use the full LUTs.
freshLithiationLimits = [min([0; negativeThetaAbs; negativeThetaRel; ...
    positiveThetaAbs; positiveThetaRel]), ...
    max([1; negativeThetaAbs; negativeThetaRel; positiveThetaAbs; positiveThetaRel])];
lithiationPlotPadding = 0.05*diff(freshLithiationLimits); % Visual margin only.
freshLithiationLimits = freshLithiationLimits+[-lithiationPlotPadding, lithiationPlotPadding];
freshSOCLimits = [-0.05, 1.05]; % Show the SOC endpoints inside the axes box.

freshOptimization = struct('positiveThetaAbsolute', positiveThetaAbs, ...
    'positiveThetaRelative', positiveThetaRel, 'positiveOCP', positiveElectrodeOCP, ...
    'negativeThetaAbsolute', negativeThetaAbs, ...
    'negativeThetaRelative', negativeThetaRel, 'negativeOCP', negativeElectrodeOCP, ...
    'absoluteEndpoints', freshThetaLimits, 'relativeEndpoints', thetaLimitsRelative(:,:,1), ...
    'lithiationLimits', freshLithiationLimits, 'socLimits', freshSOCLimits);
figures(2) = plotStateOptimization('Fresh-state optimization', cellSOC, ...
    reconstructedOCV(:,1), freshThetaLimits, thetaLimitsRelative(:,:,1), ...
    freshOptimization, gemColors);

%% Fresh and aged OCV versus discharged capacity
% SOH is usable capacity relative to the fresh state, in percent.
colorFactors = linspace(1, 0.3, nStates);
figures(3) = figure('Name', 'OCV versus capacity');
ax = axes('Parent', figures(3));
hold(ax, 'on');
for iState = 1:nStates
    plot(ax, dischargedCapacity(:,iState), reconstructedOCV(:,iState), ...
        'Color', gemColors(1,:)*colorFactors(iState), 'LineWidth', 2, ...
        'DisplayName', sprintf('%s: SOH = %.1f\\%%', ...
        stateNames(iState), capacityRetention(iState)));
end
xlabel(ax, 'Capacity (Ah)');
ylabel(ax, 'Voltage (V)');
xlim(ax, [0 max(cellCapacity)]);
legend(ax, 'Location', 'best');
grid(ax, 'on');
thesisFormat([], [], figures(3));

%% Fresh and aged OCV versus SOC
figures(4) = figure('Name', 'OCV versus SOC');
ax = axes('Parent', figures(4));
hold(ax, 'on');
for iState = 1:nStates
    plot(ax, cellSOC, reconstructedOCV(:,iState), ...
        'Color', gemColors(1,:)*colorFactors(iState), 'LineWidth', 2, ...
        'DisplayName', sprintf('%s: SOH = %.1f\\%%', ...
        stateNames(iState), capacityRetention(iState)));
end
xlabel(ax, 'SOC $z$', 'Interpreter', 'latex');
ylabel(ax, 'Voltage (V)');
xlim(ax, [0 1]);
% legend(ax, 'Location', 'best');
grid(ax, 'on');
% Match the normalized inset sizes and positions in the graphite OCP example.
mainAxes = ax;
insetSOCLimits = [0.05, 0.15; 0.4, 0.6; 0.8, 0.9];
insetPositions = [0.2, 0.75, 0.125, 0.125; ...
                  0.45, 0.5, 0.125, 0.125; ...
                  0.75, 0.425, 0.125, 0.125];
inset_ax = gobjects(size(insetSOCLimits, 1), 1);
for iInset = 1:size(insetSOCLimits, 1)
    inset_ax(iInset) = axes('Parent', figures(4), 'Units', 'normalized', ...
        'Position', insetPositions(iInset,:), 'Box', 'on');
    % Let thesisFormat associate this inset with the main plot when resizing.
    setappdata(inset_ax(iInset), 'thesisFormatMainAxes', mainAxes);
    hold(inset_ax(iInset), 'on');
    inWindow = cellSOC >= insetSOCLimits(iInset,1) ...
        & cellSOC <= insetSOCLimits(iInset,2);
    for iState = 1:nStates
        plot(inset_ax(iInset), cellSOC(inWindow), reconstructedOCV(inWindow,iState), ...
            'Color', gemColors(1,:)*colorFactors(iState), 'LineWidth', 2);
    end
    xlim(inset_ax(iInset), insetSOCLimits(iInset,:));
end
thesisFormat([], [], figures(4));
% Keep voltage tick labels legible in the fixed-size zoom axes.
for iInset = 1:numel(inset_ax)
    insetTicks = inset_ax(iInset).YTick;
    if numel(insetTicks) > 2
        inset_ax(iInset).YTick = insetTicks([1, end]);
    end
end

%% Last aged state: use the same layout and scales as the fresh optimization
figures(5) = plotStateOptimization( ...
    sprintf('Last aged-state optimization (%s)', stateNames(end)), cellSOC, ...
    reconstructedOCV(:,end), thetaLimits(:,:,end), thetaLimitsRelative(:,:,end), ...
    freshOptimization, gemColors);
coordinateAnalysis = struct('calibration', coordinateCalibrationSummary, ...
    'thetaLimitsRelative', thetaLimitsRelative, ...
    'relativeElectrodeCapacity', relativeElectrodeCapacity, ...
    'inventoryAbsoluteEndpoints', inventoryAbsoluteEndpoints, ...
    'inventoryRelativeEndpoints', inventoryRelativeEndpoints, ...
    'inventoryCorrectedEndpoints', inventoryCorrectedEndpoints, ...
    'inventoryOrigin', inventoryOrigin, 'lliComparison', lliComparison, ...
    'lamComparison', lamComparison);
% Wrap cell-valued fields in scalar cells to keep results a scalar structure.
results = struct('stateNames', {stateNames}, 'degradationModes', degradationModes, ...
    'cellSOC', cellSOC, 'thetaGuesses', thetaGuesses, ...
    'freshThetaLimits', freshThetaLimits, 'thetaLimits', thetaLimits, ...
    'freshOptimization', freshOptimization, ...
    'ocv', reconstructedOCV, 'cellCapacity', cellCapacity, ...
    'electrodeCapacity', electrodeCapacity, 'residuals', {residuals}, ...
    'exitflag', exitflag, 'summary', stateSummary, 'figures', figures, ...
    'coordinateAnalysis', coordinateAnalysis, 'paperTable', paperTable, ...
    'randomSeed', randomSeed, 'randomLossRange', randomLossRange, ...
    'randomLossStep', randomLossStep, 'randomLossLevels', randomLossLevels, ...
    'solverOptions', options, ...
    'cutOffVoltage', cutOffVoltage, 'voltageTolerance', voltageTolerance, ...
    'capacityTolerance', capacityTolerance, ...
    'extrapolation', extrapolationSummary, ...
    'referenceProvenance', referenceProvenance, ...
    'inset_ax', inset_ax);

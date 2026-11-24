function [thetaLimits, ...
          reconstructedOCV, ...
          cellCapacity, ...
          electrodeCapacity, ...
          exitflag, residual] = buildOCV(negativeElectrodeOCPfun, positiveElectrodeOCPfun, cutOffVoltage, ...
              negativeElectrodeCapacity, positiveElectrodeCapacity, thetaPrev, cellSOC, LLI, LAMn, LAMp, options)
%BUILDOCV Generate an OCV state from prescribed fresh-relative losses.
%   [thetaLimits, reconstructedOCV, cellCapacity, electrodeCapacity, ...
%       exitflag, residual] = buildOCV(Un, Up, cutOffVoltage, Qn, Qp, ...
%       thetaPrev, cellSOC, LLI, LAMn, LAMp, options).
%   Un and Up are electrode OCP interpolants; Qn and Qp are fresh capacities.
%   cutOffVoltage = [minimum maximum] [V]; options is optional fmincon options.
%   thetaLimits is 2-by-2 with the same layout as thetaPrev; cellCapacity is
%   the usable cell capacity [Ah]; electrodeCapacity = [Qn_aged Qp_aged] [Ah].
%   exitflag is the fmincon termination flag (positive indicates convergence).
%   OCP functions take absolute stoichiometry [-] and return potential [V].
%   Capacities are full-stoichiometry capacities [Ah]; losses are fractions.
%   thetaPrev is 2-by-2: rows SOC0/SOC1, columns negative/positive. For every
%   aged call it must be the SAME fitted fresh reference, not a prior age.
%   cellSOC is in [0,1]: 0 = end of discharge, 1 = end of charge.
%   reconstructedOCV is a row vector [V] in the supplied cellSOC order.
%   residual reports cutoff errors [V], balance errors [Ah], and solver data.

if nargin < 11 || isempty(options)
    options = optimoptions('fmincon', 'Display', 'off', ...
        'OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-10);
end
validateattributes(cutOffVoltage, {'numeric'}, {'vector', 'numel', 2, 'real', 'finite', 'increasing'});
validateattributes(thetaPrev, {'numeric'}, {'size', [2 2], 'real', 'finite', '>=', 0, '<=', 1});
validateattributes(cellSOC, {'numeric'}, {'vector', 'real', 'finite', '>=', 0, '<=', 1});
validateattributes(negativeElectrodeCapacity, {'numeric'}, {'scalar', 'real', 'positive', 'finite'});
validateattributes(positiveElectrodeCapacity, {'numeric'}, {'scalar', 'real', 'positive', 'finite'});
validateattributes(LLI, {'numeric'}, {'scalar', 'real', '>=', 0, '<', 1, 'finite'});
validateattributes(LAMn, {'numeric'}, {'scalar', 'real', '>=', 0, '<', 1, 'finite'});
validateattributes(LAMp, {'numeric'}, {'scalar', 'real', '>=', 0, '<', 1, 'finite'});
if thetaPrev(2,1) <= thetaPrev(1,1) || thetaPrev(2,2) >= thetaPrev(1,2)
    error('buildOCV:InvalidReference', 'Reference endpoints must have positive charge spans.');
end

cutOffVoltage = cutOffVoltage(:);
cellSOC = cellSOC(:).';
Qna = negativeElectrodeCapacity * (1 - LAMn);
Qpa = positiveElectrodeCapacity * (1 - LAMp);
% Define fresh inventory at SOC 0; LAM changes capacities, while LLI sets
% the target inventory independently relative to this same fresh reference.
Lif = negativeElectrodeCapacity * thetaPrev(1, 1) + positiveElectrodeCapacity * thetaPrev(1, 2);
Lia_Target = (1 - LLI) * Lif;

% Variables: [theta_n,EoD; theta_p,EoD; theta_n,EoC; theta_p,EoC].
xStart = [thetaPrev(1, 1), thetaPrev(1, 2), ...
          thetaPrev(2, 1), thetaPrev(2, 2)];
lb = zeros(1, 4);
ub = ones(1, 4);
% fmincon enforces Aeq*x = beq: the first row fixes inventory at EoD;
% the second keeps inventory equal at EoC, enforcing equal charge throughput.
Aeq = [ Qna,  Qpa,   0,   0;
       -Qna, -Qpa, Qna, Qpa];
beq = [Lia_Target; 0];
% Negative stoichiometry rises and positive stoichiometry falls on charge.
A = [1,  0, -1, 0;
     0, -1,  0, 1];
b = [0; 0];
% Normalize voltage residuals to mV so small squared errors do not trigger
% premature solver termination. This constant scaling preserves the optimum.
voltageScale = 1e-3; % V; numerical scale, not a measurement uncertainty.
objective = @(x) endpointResidual(x, negativeElectrodeOCPfun, ...
    positiveElectrodeOCPfun, cutOffVoltage)/voltageScale^2;
[x, ~, exitflag, output] = fmincon(objective, xStart, ...
    A, b, Aeq, beq, lb, ub, [], options);

thetaLimits = [x(1), x(2); x(3), x(4)];
QnThroughput = Qna * (x(3) - x(1));
QpThroughput = Qpa * (x(2) - x(4));
% The equality constraint makes these equal up to the solver tolerance.
cellCapacity = mean([QnThroughput, QpThroughput]);
electrodeCapacity = [Qna, Qpa];

% With constant electrode capacities, stoichiometry varies linearly with SOC.
thetaNegative = x(1) + cellSOC * (x(3) - x(1));
thetaPositive = x(2) + cellSOC * (x(4) - x(2));
reconstructedOCV = positiveElectrodeOCPfun(thetaPositive) - negativeElectrodeOCPfun(thetaNegative);

[~, voltageError] = endpointResidual(x, negativeElectrodeOCPfun, ...
    positiveElectrodeOCPfun, cutOffVoltage);
residual.voltageEoD = voltageError(1);
residual.voltageEoC = voltageError(2);
residual.capacityMismatch = QnThroughput - QpThroughput;
residual.lithiumInventory = Qna*x(1) + Qpa*x(2) - Lia_Target;
residual.voltageScale = voltageScale;
residual.optimization = output;
end

function [value, voltageError] = endpointResidual(x, Uocpn, Uocpp, Vcut)
%ENDPOINTRESIDUAL Evaluate squared cutoff errors at both SOC endpoints.
%   [value, voltageError] = endpointResidual(x, Uocpn, Uocpp, Vcut).
%   x = [theta_n0 theta_p0 theta_n1 theta_p1], dimensionless; OCPs return V.
%   Vcut = [minimum maximum] [V]. Outputs are a scalar [V^2] and 2-by-1 [V].
%   Query both endpoints together to avoid repeated interpolation preparation.
negativePotential = Uocpn(x([1 3]));
positivePotential = Uocpp(x([2 4]));
voltageError = positivePotential(:) - negativePotential(:) - Vcut(:);
value = sum(voltageError.^2);
end

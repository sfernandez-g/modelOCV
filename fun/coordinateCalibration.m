function [origin, scale] = coordinateCalibration(thetaAbsolute, thetaRelative)
%COORDINATECALIBRATION Recover a fixed affine map from paired LUT samples.
%   [origin, scale] = coordinateCalibration(thetaAbsolute, thetaRelative).
%   Inputs are equally sized, finite dimensionless vectors. Outputs satisfy
%   thetaAbsolute = origin + scale*thetaRelative, with a positive scale.
validateattributes(thetaAbsolute, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
validateattributes(thetaRelative, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
if numel(thetaAbsolute) ~= numel(thetaRelative)
    error('modelOCV:CoordinateSizeMismatch', 'Absolute and relative LUTs must contain paired samples.');
end
thetaAbsolute = thetaAbsolute(:);
thetaRelative = thetaRelative(:);
[relativeMinimum, iMinimum] = min(thetaRelative);
[relativeMaximum, iMaximum] = max(thetaRelative);
if relativeMaximum <= relativeMinimum
    error('modelOCV:DegenerateCalibration', 'Relative LUT must span a nonzero interval.');
end
% Use paired samples at the relative extrema; all pairs are checked below
% so nonlinear coordinate mappings are rejected rather than approximated.
scale = (thetaAbsolute(iMaximum)-thetaAbsolute(iMinimum))/(relativeMaximum-relativeMinimum);
origin = thetaAbsolute(iMinimum)-scale*relativeMinimum;
% Allow roundoff only: PDF Eqs. 21-26 require a genuinely affine calibration.
mappingTolerance = 1e-10*max(1, max(abs(thetaAbsolute)));
if ~isfinite(scale) || scale <= 0 || ~isfinite(origin) ...
        || max(abs(thetaAbsolute-(origin+scale*thetaRelative))) > mappingTolerance
    error('modelOCV:NonAffineCalibration', ...
        'Paired LUTs must satisfy thetaAbs = a + b*thetaRel with fixed b > 0.');
end
end

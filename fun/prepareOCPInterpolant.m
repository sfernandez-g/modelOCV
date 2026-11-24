function potential = prepareOCPInterpolant(theta, voltage)
%PREPAREOCPINTERPOLANT Prepare a linear electrode potential lookup once.
%   potential = prepareOCPInterpolant(theta, voltage); U = potential(thetaQuery).
%   theta is dimensionless and voltage is in V; both are matching finite
%   real vectors. Duplicate theta samples are averaged before interpolation.
%   Returns a griddedInterpolant with linear interpolation/extrapolation.
validateattributes(theta, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
validateattributes(voltage, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
if numel(theta) ~= numel(voltage)
    error('prepareOCPInterpolant:SizeMismatch', 'theta and voltage must have matching sizes.');
end
% Apply the same permutation to both vectors to preserve sample pairing.
[theta, order] = sort(theta(:));
voltage = voltage(order);
[theta, ~, group] = unique(theta);
if numel(theta) < 2
    error('prepareOCPInterpolant:InsufficientData', 'At least two distinct theta values are required.');
end
% unique assigns a group index to each sample; accumarray averages its voltages.
voltage = accumarray(group, voltage(:), [], @mean);
% The two methods specify interpolation inside and extrapolation outside the LUT.
potential = griddedInterpolant(theta, voltage, 'linear', 'linear');
end

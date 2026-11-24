# modelOCV

MATLAB example accompanying a scientific article on lithium-inventory loss
(LLI) and active-material loss (LAM) under absolute and relative electrode
lithiation coordinates. The synthetic scenarios illustrate coordinate
consistency; they do not estimate degradation from experimental aged-cell data.

## Requirements and use

MATLAB R2023b or newer and Optimization Toolbox.

```matlab
run('modelOCV.m')
```

The script locates `fun/` and `input/` relative to itself. It clears its execution
workspace, closes existing figures, displays tables, opens five figures and
leaves `results` and intermediate variables in the workspace. It saves no files.
From another directory, use `run('/path/to/modelOCV/modelOCV.m')`.

## Files

- `modelOCV.m`: model parameters, synthetic scenarios, calculations and figures.
- `fun/`: OCV construction, OCP interpolation, coordinate calibration and plotting.
- `input/Uocp_n.mat`: negative-electrode reference model.
- `input/Uocp_p.mat`: positive-electrode reference model.

The local `doc/` and independent `modelOCP/` directories are ignored by Git and
are not required to run this example.

## Reference data

The script reads only the numerical average-branch arrays `model.LUT.theta_avg`
(absolute stoichiometry), `thetarel_avg` (relative lithiation) and `U_avg`
(potential versus Li, V) from each supplied MAT file. No fitting, resampling or
changes to the source files are performed. The files also contain legacy model
data and serialized function handles; those handles are not used. MATLAB may
warn about unavailable functions when loading them.

The source SHA-256 identifiers are:

| File | SHA-256 |
| --- | --- |
| Uocp_n.mat | B6C5ADD6D00953E0213AE7AA6FDC2B3FED13C896CF572CA74E733426CCF8E885 |
| Uocp_p.mat | 2FD8919F4DED206B4C7A5762E00B465AD4498B8B969419F0A6868E82D14E801D |

These supplied references define the article's calibration. Regenerated outputs
from the companion [modelOCP](https://github.com/sfernandez-g/modelOCP) project
are not automatically interchangeable with them. Absolute stoichiometry remains
conditional on the original electrode-calibration assumptions.

## Numerical conventions

- SOC 0 is end of discharge and SOC 1 is end of charge. Endpoint matrices have
  rows `[SOC0; SOC1]` and columns `[negative positive]`.
- Capacities and lithium inventories are charge equivalents in Ah. Losses are
  fractions internally and percentages in the displayed tables.
- Every aged state uses the same fresh capacities and fitted fresh inventory.
- The fixed affine map is `thetaAbsolute = origin + scale*thetaRelative`;
  relative electrode capacities equal `scale*absoluteCapacity`.
- OCP interpolation and extrapolation are linear. `results.extrapolation`
  reports excursions outside the sampled domains.
- A local random stream with seed 2027 defines the mixed degradation scenarios.
  Inputs, solver settings and residuals are retained in `results`.

The script checks solver convergence, finite values, positive operating spans,
cutoff residuals and charge/inventory balances. These checks do not establish
experimental accuracy or the physical validity of extrapolation.

## License and citation

See [LICENSE.md](LICENSE.md) for the code license (CC BY-SA 4.0), numerical LUT
license (CC BY 4.0), and scope exclusions. When citing this software, identify
the source revision and input data used. No release DOI is assigned here.

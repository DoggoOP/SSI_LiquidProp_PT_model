# SSI Liquid Propellant Rocket Engine - Regenerative Cooling Analysis

A MATLAB-based thermal and fluid analysis tool for liquid rocket engines with regenerative cooling. This project performs detailed heat transfer calculations, pressure drop analysis, and nozzle geometry generation for rocket engine design.

## Overview

The overarching goal of this project is to create a Nitrous IPA liquid propellant rocket. This project models the regenerative cooling system of a liquid rocket engine, where the fuel flows through cooling channels around the combustion chamber and nozzle before being injected. The analysis calculates:

- Temperature distribution along the engine (gas, wall, and coolant temperatures)
- Heat transfer coefficients (gas-side and coolant-side)
- Pressure drop through cooling channels
- Required feed system pressure
- Thermal and hydraulic performance metrics

## Files

### RegenCooling.m
Main analysis script that performs:
- Nozzle geometry generation (cylindrical chamber, converging, and diverging sections)
- Isentropic flow calculations (Mach number, area ratio)
- Iterative heat transfer analysis using:
  - Bartz equation for gas-side heat transfer
  - Gnielinski correlation for coolant-side heat transfer
  - Fin efficiency corrections for channel geometry
- Pressure drop calculations using Darcy-Weisbach equation
- Comprehensive plotting and results output

**Key Parameters:**
- Coolant: Fuel (likely RP-1 based on properties)
- Chamber pressure: 2.068 MPa
- Stagnation temperature: 2530 K
- Throat diameter: 24.36 mm
- Exit diameter: 44.92 mm
- Mass flow rate: 0.168 kg/s
- Number of cooling channels: 35
- Channel dimensions: 1.5 mm x 1.5 mm

### NozzleContour.m
Generates bell-shaped nozzle contour using the Rao method.

**Function Signature:**
```matlab
r = NozzleContour(dt, de, theta_n, theta_e, points)
```

**Parameters:**
- `dt` - Throat diameter (m)
- `de` - Exit diameter (m)
- `theta_n` - Nozzle wall match angle (degrees)
- `theta_e` - Exit angle (degrees)
- `points` - Number of discretization points

**Returns:**
- `r` - Radius as a function of axial position

The contour consists of:
1. Circular arc near the throat (radius of curvature = 1.5 × throat radius)
2. Parabolic expansion section

### colebrook.m
Calculates the Darcy-Weisbach friction factor for pipe flow using the Colebrook-White equation.

**Function Signature:**
```matlab
f = colebrook(Re, epsilon)
```

**Parameters:**
- `Re` - Reynolds number (scalar or array)
- `epsilon` - Relative roughness (ε/D)

**Returns:**
- `f` - Friction factor

Supports vectorized operations for multiple Reynolds numbers and/or roughness values.

## Features

### Heat Transfer Analysis
- Gas-side convection using Bartz equation with compressibility corrections
- Coolant-side convection using Gnielinski correlation for turbulent flow
- Temperature-dependent coolant properties (density, viscosity, thermal conductivity, specific heat)
- Fin efficiency corrections for rectangular cooling channels
- Adiabatic wall temperature calculation with Prandtl number recovery factor

### Pressure Drop Analysis
- Darcy-Weisbach equation for friction losses
- Reynolds number-dependent friction factor:
  - Laminar flow (Re < 2300): f = 64/Re
  - Turbulent flow (Re ≥ 2300): Haaland approximation
- Surface roughness consideration (2 μm for copper)
- Required feed pressure calculation based on injector pressure requirements

### Flow Analysis
- Isentropic flow relations for entire chamber and nozzle
- Mach number distribution from chamber (M ≈ 0) to supersonic exit
- Area ratio calculations
- Mass flow distribution across cooling channels

## Dependencies

### Required MATLAB Toolboxes:
- Aerospace Toolbox (for `flowisentropic` function)
- Optimization Toolbox (for `fzero` function in colebrook.m)

### Standard Functions Used:
- Plotting: `plot`, `surf`, `figure`, `legend`, etc.
- Mathematical: `sqrt`, `exp`, `log10`, `tanh`
- Array operations: `linspace`, `meshgrid`, `reshape`

## Usage

### Running the Main Analysis

```matlab
% Simply run the main script
RegenCooling
```

The script will:
1. Generate the nozzle contour
2. Calculate flow properties throughout the engine
3. Perform iterative heat transfer calculations
4. Compute pressure drop
5. Display results in console
6. Generate multiple analysis plots

### Generating a Custom Nozzle Contour

```matlab
% Parameters
dt = 24.36E-3;      % Throat diameter (m)
de = 44.92E-3;      % Exit diameter (m)
theta_n = 17.96;    % Match angle (degrees)
theta_e = 8;        % Exit angle (degrees)
points = 100;       % Number of points

% Generate contour
r = NozzleContour(dt, de, theta_n, theta_e, points);

% Plot
x = linspace(0, length_nozzle, points);
plot(x, r)
```

### Calculating Friction Factor

```matlab
% Single value
Re = 1e5;
epsilon = 1e-4;
f = colebrook(Re, epsilon);

% Multiple Reynolds numbers
Re_range = logspace(4, 8, 100);
f_values = colebrook(Re_range, epsilon);
plot(Re_range, f_values)
```

## Output

### Console Output
The script prints:
- Total pressure drop (in MPa, Pa, psi, bar)
- Required injector pressure
- Required feed pressure
- Reynolds number statistics (mean, min, max)
- Friction factor statistics
- Channel geometry parameters
- Maximum wall temperature and location
- Design summary

### Generated Plots
1. Nozzle radius contour
2. Complete chamber/nozzle geometry
3. Gas-side heat transfer coefficient distribution
4. Coolant heat transfer coefficient distribution
5. Temperature distribution (wall, coolant, contour)
6. Adiabatic wall temperature
7. Heat flux distribution
8. Coolant pressure distribution
9. Local pressure drop per segment
10. Dual-axis pressure and velocity plot
11. Reynolds number distribution
12. Friction factor distribution

## Methodology

### Heat Transfer
- **Bartz Equation**: Industry-standard correlation for gas-side heat transfer in rocket nozzles, accounting for Mach number effects, throat curvature, and area ratio
- **Gnielinski Correlation**: Accurate for turbulent flow in channels with moderate Prandtl numbers
- **Fin Theory**: Cooling channels approximated as rectangular fins with efficiency corrections

### Pressure Drop
- **Darcy-Weisbach**: Fundamental equation for pressure drop in pipes
- **Colebrook-White**: Implicit equation for friction factor in turbulent flow with surface roughness
- **Haaland Approximation**: Explicit approximation to Colebrook equation

### Flow Modeling
- **Isentropic Relations**: Assumes ideal gas with constant gamma, valid for preliminary design
- **Subsonic/Supersonic Regions**: Properly handled using MATLAB's `flowisentropic` function

## Design Notes

### Assumptions
1. Isentropic flow throughout chamber and nozzle
2. Constant specific heat ratio (γ = 1.24)
3. Steady-state operation
4. Negligible heat conduction in axial direction
5. Uniform flow distribution across all cooling channels
6. Temperature-dependent coolant properties (linearized approximations)

### Convergence Criteria
- Wall temperature iteration converges when ΔT < 15 K between iterations
- Ensures accurate coupling between gas-side and coolant-side heat transfer

### Design Variables to Modify
Key parameters that can be adjusted in `RegenCooling.m`:
- Line 6: `T1` - Coolant inlet temperature
- Line 9: `mf` - Fuel mass flow rate
- Line 19: `Pc` - Chamber pressure
- Line 22-23: `d`, `w` - Channel dimensions
- Line 110: `N` - Number of cooling channels
- Line 134: `P_injector_required` - Required injector pressure

## References

### Textbooks/Literature
- Huzel, D. K., & Huang, D. H. (1992). *Modern Engineering for Design of Liquid-Propellant Rocket Engines*
- Sutton, G. P., & Biblarz, O. (2016). *Rocket Propulsion Elements*

### Correlations
- Bartz, D. R. (1957). "A Simple Equation for Rapid Estimation of Rocket Nozzle Convective Heat Transfer Coefficients"
- Colebrook, C. F. (1939). "Turbulent Flow in Pipes, with Particular Reference to the Transition Region"
- Gnielinski, V. (1976). "New Equations for Heat and Mass Transfer in Turbulent Pipe and Channel Flow"

## Author

SSI (Students for the Exploration and Development of Space)
Liquid Propellant Rocket Engine Team

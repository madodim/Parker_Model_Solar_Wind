# Parker_Model_Solar_Wind
A MATLAB-based model of the Parker solar wind solution with physical unit conversion, density profiles, and comparison to empirical coronal models.

**Solar Wind Model - Parker Equation and Solar Corona Dynamics**\\
This repository contains MATLAB scripts for modeling the solar wind based on the Parker solar wind equation. It includes visualization of isodense contours, selection of physically relevant wind solutions, interpolation, transformation into real physical units (velocity, distance, density), and comparison with empirical solar atmospheric models.

**Code Overview**\\
1. Isodense Contours\\
Plots contours of constant "C" values from the dimensionless Parker solar wind equation:
- Shows the structure of possible solutions in the (ξ, u) space.
- Marks the critical point where the velocity equals the speed of sound.

2. Solution Choice\\
Extracts and filters the physically valid solution for a chosen isodense value (C = -3), selecting points that start subsonic and end supersonic.

3. Interpolation\\
Interpolates the selected Parker wind solution using linear interpolation to create a smooth solar wind velocity profile in normalized units.

4. Velocity and Distance (Physical Units)
- Converts the normalized Parker wind solution into physical units (km/s and solar radii).
- Plots solar wind speed for various coronal temperatures.
- Overlays Earth’s orbit for context.

5. Density vs Distance
- Calculates and plots the solar wind plasma density in cm⁻³ as a function of distance.
- Uses the mass conservation law and outputs results for multiple coronal temperatures.
- Logarithmic scale captures the density drop-off.

6. Save Density and Frequency Profiles
- Saves the density and plasma frequency (Hz) at selected heliocentric distances (in solar radii) to .txt files. Useful for radio burst and solar wind diagnostics.

7. Atmospheric Models\\
Implements widely used empirical models: Newkirk, Saito, Leblanc, and Vršnak. Compares their predicted plasma densities across heliocentric distances.

**Files**
- parker_solar_wind.m: Main MATLAB script containing all sections described above.
- Output files with plasma density and frequency profiles. Change the file name and the directory accordingly.
- README.md: Documentation file (you’re reading it!).

**📌 Requirements**\\
MATLAB (R2018+ recommended)
No external toolboxes required

**📘 References**
- Parker, E. N. (1958). Dynamics of the interplanetary gas and magnetic fields. Astrophysical Journal, 128, 664.
- Newkirk, G. Jr. (1961). The Solar Corona in Active Regions and the Thermal Origin of the Slowly Varying Component of Solar Radio Radiation. Astrophysical Journal, 133, 983.
- Saito, K., Poland, A. I., & Munro, R. H. (1977). A Study of the Background Corona Near Solar Minimum. Solar Physics, 55, 121–134.
- Leblanc, Y., Dulk, G. A., & Bougeret, J.-L. (1998). Tracing the Electron Density from the Corona to 1 AU. Solar Physics, 183, 165–180.
- Vršnak, B. et al. (2004). Band-splitting of coronal and interplanetary type II bursts. III. Physical conditions in the upper corona and interplanetary space. A&A, 413, 753–763.

% Ekin Gunes Ozaktas, June 2024
%
% Code to plot extinction profile for periodic square lattice of plasmonic
% nanodisks.
%
% For more details see E. G. Ozaktas et al. "Aperiodicity and Disorder as
% Systematic Spectral Tuning Mechanisms for Plasmonic Nanostructures",
% 2025.

clear all;

h = 6.626e-34; % Planck constant
e = 1.602e-19; % Elementary charge
c = 2.998e8; % Speed of Light

% Data file
load("BasicPhysics/BasicPhysics_a_1200nm_r100_h200.mat");

FS = 20; % Font Size

% Plotting

figure()
plot(h*f/e, -log(T), "LineWidth", 1.5);
xlabel("Energy (eV)")
ylabel("Extinction")
fontsize(FS,"points")
axis([min(h*f/e), max(h*f/e), 0, 0.5])

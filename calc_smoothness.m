% Ekin Gunes Ozaktas, June 2024
%
% Calculation of local variance for a given frequency spectrum, with width
% Df.
%
% For more details see E. G. Ozaktas et al. "Aperiodicity and Disorder as
% Systematic Spectral Tuning Mechanisms for Plasmonic Nanostructures",
% 2025.

function smoothness = calc_smoothness(f, data, Df)

    % Convert frequency width to index width
    df = f(2)-f(1);
    N = round(Df/df);
    
    % Moving average
    smoothed = smooth(data, N);
    
    var = (data - smoothed).^2 ./ smoothed.^2;
    
    % Integrate
    smoothness = sum(var)*df;
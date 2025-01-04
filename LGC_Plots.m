% Ekin Gunes Ozaktas, June 2024
%
% Code to plot extinction profile, calculate local variance, and visualize geometry 
% for plasmonic nanodisks distributed according to a log Gaussian Cox process.
%
% For more details see E. G. Ozaktas et al. "Aperiodicity and Disorder as
% Systematic Spectral Tuning Mechanisms for Plasmonic Nanostructures",
% 2025.

clear all;
h = 6.626e-34; % Planck constant
e = 1.602e-19; % Elementary charge
c = 2.998e8; % Speed of Light

% Calculation of parameter values and loading file names
NUMTRIALS = 5;
params = [0.5,1,1.5,2,2.5,3,4,5];
param_names = ["05", "1", "15", "2", "25", "3", "4", "5"];
P_VEC = params;
L = 4.8e-6;
grid_size = [30, 30];
a_vec = L/grid_size(1);
fnamesbase = "LGC/LGC_s" + param_names + "_l02_N30";
fnames_list = fnamesbase + ".mat";
for ind = 1:length(fnames_list)
    names(ind) = "s = " + num2str(P_VEC(ind));
end

r = 100e-9; % radius in nm
startf = 237; % frequency index to start from when taking local variance

newcolors = {'#FF0000','#FF8000','#999900', '#80FF00', '#00FF00', '#0080FF', '#0000FF', '#7F00FF'};
FS = 20; % Font size

Df = 0.7; % Moving average window width

% Plot spectra and calculate local variances.
figure()
hold on
voffset = 0.2;
for ind = 1:length(fnames_list)
load(fnames_list(ind));
smoothness_range = startf:length(f);
plot(h*f/e, -log(squeeze(T)) + (length(fnames_list) - ind)*voffset, "LineWidth", 1.5);

load(fnamesbase(ind) + ".mat");
        smoothness_range = startf:length(f);
        smoothnessE(ind,1) = calc_smoothness(h*f(smoothness_range)/e, -log(squeeze(T)), Df);    
    for ind2 = 2:(NUMTRIALS)
        load(fnamesbase(ind) + "_" + num2str(ind2) + ".mat");
        smoothness_range = startf:length(f);
        smoothnessE(ind,ind2) = calc_smoothness(h*f(smoothness_range)/e, -log(squeeze(T)), Df);
    end

end
legend(names)
colororder(newcolors)
xlabel("Energy (eV)")
ylabel("Extinction")
axis([min(h*f/e), max(h*f/e), 0, 1.6])
fontsize(FS,"points")

% Plot local variance for all spectra.
figure()
hold on
for ind = 1:length(fnames_list)
scatter(P_VEC(ind), smoothnessE(ind,1), "LineWidth", 1.5)
end
for ind2 = 2:NUMTRIALS
for ind = 1:length(fnames_list)
scatter(P_VEC(ind), smoothnessE(ind,ind2), "LineWidth", 1.5)
end
end
colororder(newcolors)
xlabel("s")
ylabel("Local Variance")
fontsize(FS,"points")

% Plot geometric distributions.

figure()
hold on
for ind2 = 1:length(fnames_list)
    
    load(fnames_list(ind2));
    
    a = 0.5e-6;
    L = 4.8e-6;
    kx = 0:(2*pi/L):(2*pi/a - 2*pi/L);
    ky = 0:(2*pi/L):(2*pi/a - 2*pi/L);
    [KX, KY] = meshgrid(kx,ky);
    FT = 0*KX;
    
    subplot(4,length(fnames_list),ind2)
    hold on
    color_gold = 'y';
    % Plot real space and calculate Fourier Transform
    for l = 1:length(locs)
        x = (L-2*r)*locs(l,1) + r;
        y = (L-2*r)*locs(l,2) + r;
        FT = FT + exp(-1j*(KX*x + KY*y));
        p = nsidedpoly(1000, 'Center', [x, y]*1e6, 'Radius', r*1e6);
        plot(p, 'FaceColor', color_gold)
    end
    title(names(ind2))
    axis equal;
    axis([0,grid_size(1)*a_vec,0,grid_size(2)*a_vec]*1e6)
    ax = gca;
    ax.FontSize = 12;
    ax.TitleFontSizeMultiplier = 1.4;
    ax.LabelFontSizeMultiplier = 1.4;
    box on
    ax.XTick = [0,L]*1e6;
    ax.XTickLabel = {"0",sprintf("%0.2f",1e6*L)};
    ax.YTick = [0,L]*1e6;
    ax.YTickLabel = {"0",sprintf("%0.2f",1e6*L)};
    
    win = 1;
    res_mat = zeros(win, win);
    res_mat(ceil(win/2), ceil(win/2)) = 1;
    FT_img = kron(FT, res_mat);
    
    subplot(4,length(fnames_list),length(fnames_list) + ind2)
    hold on
    colormap hot;
    axis equal;
    sz = size(FT_img);
    axis([1-0.5,sz(1)+0.5,1-0.5,sz(2)+0.5])
    
    % Plot Fourier Transform
    xl = 1:sz(1);
    yl = 1:sz(2);
    [X,Y] = meshgrid(xl,yl);
    X = [reshape(X, [sz(1)*sz(2),1]); 0];
    Y = [reshape(Y, [sz(1)*sz(2),1]); 0];
    F = [reshape(abs(FT_img), [sz(1)*sz(2),1]); 0];
    [sortF, sortinds] = sort(F);
    scatter(X(sortinds),Y(sortinds),[],F(sortinds),'filled','LineWidth',0.1);
    ax = gca;
    ax.Color = 'k';
    
    set(gca,'TickDir','out');
    xticks([1, sz(1)*win]);
    yticks([1, sz(2)*win]);
    xticklabels({"0", num2str(2*(sz(1)*win - 1)) + "\pi/L"});
    yticklabels({"0", num2str(2*(sz(2)*win - 1)) + "\pi/L"});
    set(gca,'TickDir','out');
    ax = gca;
    ax.TickLength = [0.05, 0.05]; 
    
    ax = gca;
    ax.FontSize = 12;
    ax.TitleFontSizeMultiplier = 1.4;
    ax.LabelFontSizeMultiplier = 1.4;

end


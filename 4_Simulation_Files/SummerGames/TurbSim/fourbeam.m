% Four-Beam Lidar Coherence Calculation for Spectral Estimation
% Simplified to produce only a coherence vs wave number plot

clear; close all; clc;

%% Configuration Parameters

% IEA 3.4 MW Turbine Parameters
HubHeight     = 110;   % Hub height (m)
RotorDiameter = 130;   % Rotor diameter (m)

% Lidar Configuration (Four-beam pulsed lidar)
focal_distance = 100;              % Focal distance (m)
half_cone_angle = 15;              % Half cone angle (degrees)
beam_positions = [0, 90, 180, 270];% Beam azimuth angles (degrees)

% Analysis Parameters
fs      = 4;    % Sampling frequency (Hz) - typical for lidar
nfft    = 256;  % FFT length for spectral estimation
overlap = 0.5;  % 50% overlap for Welch method

%% Load / Generate Data

fprintf('Loading TurbSim data...\n');
turbsim_file = 'TurbSim2a.inp';

if exist(turbsim_file, 'file')
    fid = fopen(turbsim_file, 'r');
    if fid == -1
        error('Cannot open TurbSim input file');
    end
    fclose(fid);
    fprintf('TurbSim input file found: %s\n', turbsim_file);
else
    warning('TurbSim input file not found. Using default parameters.');
end

[data_exists, U_los] = load_turbsim_output();

if ~data_exists
    % Generate synthetic data for demonstration
    fprintf('Generating synthetic lidar data...\n');
    duration  = 600;           % 10 minutes
    t         = (0:1/fs:duration-1/fs)';
    n_samples = length(t);

    % Simulate line-of-sight velocities for 4 beams
    U_mean_nom = 10;           % Nominal mean wind speed (m/s)
    U_los      = zeros(n_samples, 4);

    for i = 1:4
        noise  = randn(n_samples, 1);
        U_turb = filter(1, [1, -0.95], noise) * 1.5; % AR(1) process
        U_los(:, i) = U_mean_nom + U_turb;
    end
else
    [n_samples, ~] = size(U_los);
    t = (0:n_samples-1)'/fs;
end

[n_samples, n_beams] = size(U_los);
fprintf('Data loaded: %d samples, %d beams, %.1f seconds\n', ...
        n_samples, n_beams, n_samples/fs);

%% Beam separation (for info, not directly needed in plot)

theta = half_cone_angle * pi/180;   % Convert to radians
r_y   = 2 * focal_distance * sin(theta); % Lateral separation
r_z   = 0;                          % Vertical separation (for horizontal scans)

fprintf('Beam separation: %.2f m (lateral)\n', r_y);

%% Spectral Estimation using Welch Method

fprintf('Computing power spectra...\n');

window   = hann(nfft);
noverlap = round(nfft * overlap);

% Single-beam PSDs
S = zeros(nfft/2+1, n_beams);
for i = 1:n_beams
    [S(:,i), f] = pwelch(U_los(:,i), window, noverlap, nfft, fs);
end

%% Cross-Spectral Densities

fprintf('Computing cross-spectra...\n');

[Sxy_13, f] = cpsd(U_los(:,1), U_los(:,3), window, noverlap, nfft, fs); % opposite
[Sxy_24, ~] = cpsd(U_los(:,2), U_los(:,4), window, noverlap, nfft, fs); % opposite
[Sxy_12, ~] = cpsd(U_los(:,1), U_los(:,2), window, noverlap, nfft, fs); % adjacent
[Sxy_23, ~] = cpsd(U_los(:,2), U_los(:,3), window, noverlap, nfft, fs); % adjacent

%% Coherence Functions

fprintf('Computing coherence functions...\n');

Coh_13 = abs(Sxy_13).^2 ./ (S(:,1) .* S(:,3)); % opposite beams 1-3
Coh_24 = abs(Sxy_24).^2 ./ (S(:,2) .* S(:,4)); % opposite beams 2-4
Coh_12 = abs(Sxy_12).^2 ./ (S(:,1) .* S(:,2)); % adjacent beams 1-2
Coh_23 = abs(Sxy_23).^2 ./ (S(:,2) .* S(:,3)); % adjacent beams 2-3

% Mean wind speed from data (for k conversion)
U_mean = mean(U_los(:));

%% Coherence vs wave number plot (your requested figure)

% Convert frequency to wave number: k = 2*pi*f / U_mean
k = 2*pi*f / U_mean;

% Choose which coherence to plot (here beams 1-3; change as needed)
gamma_Sq_RL = Coh_13;

% Find wave number where coherence crosses ~0.5 (MCB)
gamma_target = 0.5;
[~, idx_half] = min(abs(gamma_Sq_RL - gamma_target));
MCB = k(idx_half);

% Coherence figure (your style)
figure('Name','Coherence');
hold all; grid on; box on;

plot(k, gamma_Sq_RL, 'b-', 'LineWidth', 1.5);
plot([1e-3 1e0], [0.5 0.5], 'k--', 'LineWidth', 1.0);
plot(MCB, 0.5, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

xlim([1e-3 1e0]);
set(gca, 'XScale', 'log');
xlabel('wave number [rad/m]');
ylabel('coherence [-]');

legend('Coherence 1-3','\gamma = 0.5','MCB','Location','best');

fprintf('Estimated MCB wave number: %.4g rad/m\n', MCB);

%% Helper function to load TurbSim output
function [success, U_los] = load_turbsim_output()
    % Placeholder loader; set success=false to use synthetic data
    success = false;
    U_los   = [];
end

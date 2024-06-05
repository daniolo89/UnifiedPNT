clc
clear all;
close all;
addpath("Functions\")
aaa = 1;
%% Define configuration parameters
c              = 3e8 ;                                                        % Speed of light (m/s)
f              = 30e9;                                                        % Carrier frequency (GHz)
lambda         = c / f;                                                       % Wavelength (m)
TxPower_dBm    = 20;                                                          % transmit power in dBm 
TxPower        = (10^(TxPower_dBm/10))*1e-3;                                  % transmit power in watt
NPSD           = 1e-3*10^(-174/10)*1e9;                                       % noise spectral density, NPSD = kT
Noise_Factor   = 10^(8/10);                                                   % noise factor
NoiseSampleNum = 1000;                                                        % Number of noise iterations
PointNum       = 30;                                                          % Number of points UE to generate

%% Geometry
Mr             = 64;                                                          % Number of RIS elements in a row
Mc             = 64;                                                          % Number of RIS elements in a column
risElementDist = 0.5 * lambda;                                                % RIS element distance
r_vec          = 10.^(linspace(0, log10(35), PointNum));                      % The vector r_vec contains the radial distances of the points from the origin, a logarithmically spaced vector ranging from 1 to 35. 
xyz            = [-r_vec/sqrt(2); r_vec/sqrt(2); -10*ones(1, length(r_vec))]; % XYZ coordinates of the UE points
uePos_prim     = [xyz(1,:).', xyz(2,:).', xyz(3,:).'];                        % UE positions
uePos          = uePos_prim(1,:);
bsPos          = [5, 5, 0];                                                   % Base station position
risPos         = [0, 0, 0];                                                   % RIS position

%% Signal generation parameters
Nfft           = 4096;                                                        % FFT Size
N_ftt_tilde    = 2^8;                                                         % IFFT length for angle estimation
Df             = 120e3;                                                       % Subcarrier bandwidth (GHz)
Nsc            = 3e3;                                                         % Number of subcarriers
Ts             = 1/Df;                                                        % Symbol duration
Ncp            = 288;                                                         % length in samples of CP
Nsym           = 256;                                                         % # of ofdm symbols
FreqOffset     = 0;                                                           % in Hz
sqEs           = sqrt(TxPower / Nsc);                                         % Square root of symbol energy
sqN            = sqrt(Noise_Factor * NPSD * Df / 2);                          % Square root of noise power (P_N = k T B F, where k is Boltzmann's constant, T is the standard temperature in Kelvin, B is the bandwidth, F is the noise figure. 
snr            = TxPower / (Nsc * Noise_Factor * NPSD * Df);                  % SNR in linear scale
Dt             = 0; %rand*(1/Df);                                             % Random time delay, uniformly drawn from the interval [0, 1/∆f ). 
Delta_t        = 1/(Nfft * Df);                                               % Artificial delay is constrained by the range [0, 1/(Nfft Df)].

%=======================================================================================================%
% Calculate distances
BsUeDist       = sqrt(sum((bsPos - uePos).^2, 2)); % Distance from BS to UE
BsRisDist      = sqrt(sum((bsPos - risPos).^2, 2)); % Distance from BS to RIS
RisUeDist      = sqrt(sum((risPos - uePos).^2, 2)); % Distance from RIS to UE

%=======================================================================================================%

%% Compute LoS signal 
% LoS path delay
tau_b         = (BsUeDist/c) + Dt;

% Compute channel gain gb,based on Friis' formula (using unit directivity for BS, UE). 
% Their phases are set randomly between [0,2π)
% Friis' transmission formula for the channel gain is: g=sqrt((λ/(4πd))2^)
% Compute the amplitude
g_b_amplitude = sqrt((lambda / (4 * pi * BsUeDist))^2);  
% Set random phases
g_b_phase     = 2 * pi * rand(size(BsUeDist));  % Random phase in [0, 2*pi)
% Combine amplitude and phase to get the complex channel gain
gb            = g_b_amplitude * exp(1j * g_b_phase); 

% Compute frequency offset terms
db            = exp(-1j * 2 * pi * Df * tau_b  * (0:Nsc-1));

% LoS signal 
Y_b           = gb * sqEs * db.' * ones(1, Nsym) ;

%=======================================================================================================%

%% Compute reflected signal
% Reflected path delay
tau_r      = (RisUeDist+BsRisDist)/c + Dt;

% Compute the phase profile matrix
% Generate random phase shifts for RIS elements over time
% Each element of the matrix represents a random value uniformly distributed between 0 and 1. 
% Each column represents a different time instance.
% The phase angles are expressed as complex values on the unit circle.
ris_phases = NaN (Mr*Mc,Nsym);
b_k        = 2*pi*rand(Mr*Mc,Nsym/2);
for k = 0: (Nsym/2)-1
    if k == 0
        ris_phases(:,1)   = b_k(:,1);
        ris_phases(:,2)   = -b_k(:,1);
    else
        ris_phases(:,2*k+1)   = b_k(:,k+1);
        ris_phases(:,2*k+2)   = -b_k(:,k+1);
    end
end

% Compute azimuth angle (angle in x-y plane)
az_theta   = (atan2((bsPos(2) - risPos(2)), (bsPos(1) - risPos(1)))); % AOA in azimuth from the BS to the RIS
az_phi     = (atan2((risPos(2) - uePos(2)), (risPos(1) - uePos(1)))); % AOD in azimuth from the RIS to the UE

% Compute elevation angle (angle from x-y plane to direction vector)
el_theta   = acos((bsPos(3) - risPos(3))/BsRisDist); % AOA in elevation from the BS to the RIS
el_phi     = acos((risPos(3) - uePos(3))/RisUeDist); % AOD in elevation from the RIS to the UE

% Calculate wavenumber vector
k_theta     = cal_Kv(sin(el_theta), cos(az_theta), sin(az_theta), cos(el_theta), lambda);
k_phi       = cal_Kv(sin(el_phi), cos(az_phi), sin(az_phi), cos(el_phi), lambda);


% Compute RIS element's location
risElementLoc = computeRISPositions(Mr, Mc, risElementDist);

% Calculate ur (steering vector) for the reflected path
a_theta       = exp(1j*k_theta*risElementLoc.');
a_phi         = exp(1j*k_phi*risElementLoc.');

ur            = a_theta.*a_phi *ris_phases;

% Compute channel gains  gr based on Friis' formula (using unit directivity for BS, UE and RIS elements). 
% Their phases are set randomly between [0,2π)
% Friis' transmission formula for the channel gain is: g=sqrt((λ/(4πd))2^)
% Compute the amplitude
g_r_amplitude  = sqrt((lambda / (4 * pi * BsRisDist))^2 * (lambda / (4 * pi * RisUeDist))^2);  % Reflected path gain amplitude
% Set random phases
g_r_phase      = 2 * pi * rand(size(BsUeDist));  % Random phase in [0, 2*pi)
% Combine amplitude and phase to get the complex channel gains
gr             = g_r_amplitude * exp(1j * g_r_phase);  

% Compute frequency offset terms
dr             = exp(-1j * 2 * pi * Df * tau_r  * (0:Nsc-1));

% Reflected signal 
Y_r            = gr * sqEs * ((dr.' * ones(1, Nsym)) .* (ones(Nsc,1) * ur));

%=======================================================================================================%

%% Compute received signal
% Received signal without noise
Y_NL = Y_b+Y_r;

% Received signal with noise
Y = Y_NL ; %+ sqN * (randn(size(Y_NL)) + 1j * randn(size(Y_NL)));

%=======================================================================================================%

%% Estimation of tau_b
yc            = sum(Y, 2);

% Calculate the IFFT matrix F
F             = exp(2j*pi*(0:Nfft-1)' * (0:Nsc-1) / Nfft) / Nfft;

% Set delta = 0 and calculate yc_bar
delta         = 0;
phase_shift   = exp(-1j * 2 * pi * delta * (0:Nsc-1)' * Df);
yc_delta      = yc .* phase_shift;
y_c_bar       = F * yc_delta;

% Find the index k that maximizes the absolute value
[y0, k_tilde] = max(abs(y_c_bar));

%Refine the estimate of delta, where e_k is vector comprising all zeros, except with a 1 in the kth entry estimated previously. 
e_k           = zeros(Nfft, 1);
e_k(k_tilde)  = 1;

% Define the function to add artificial delay δ 
y_hat_fun     = @(delta) -abs(e_k' * (F * (yc .* exp(-1j * 2 * pi * delta * (0:Nsc-1)' * Df)))) / y0;

% Optimization options for the quasi-Newton method using fminunc
options       = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10);
delta_tilde   = fminunc(y_hat_fun, 0, options);

% Estimate tau_b
tau_b_est_coarse  = (k_tilde - 1) * Delta_t;
tau_b_est_fine    = tau_b_est_coarse - delta_tilde;

%=======================================================================================================%

%% Estimate of tau_r
% Calculate the BS contribution d_tau_hat
d_tau_hat       = exp(-1j * 2 * pi * Df * (-tau_b_est_fine) * (0:Nsc-1)).';

% Compute gb_hat
gb_hat          = sum(yc .* d_tau_hat) / (Nsc * Nsym);

% Remove the BS contribution from the received signal Y
Yr_hat          = Y - gb_hat  .* (exp(-1j * 2 * pi * Df * (tau_b_est_fine) * (0:Nsc-1)).' * ones(1,Nsym)) ;

Yr_tilde        = sum(Yr_hat,2);

% Perform  the IFFT of Yr_tilde
y_r_bar         = F * Yr_tilde;

% Find the index k_hat that maximizes the norm of the k-th row of y_r_bar
[y0_hat, k_hat] =max(vecnorm(y_r_bar,2,2));

%Refine the estimate of delta, where e_k_hat is vector comprising all zeros, except with a 1 in the kth entry estimated previously. 
e_k_hat         = zeros(Nfft, 1);
e_k_hat(k_hat)  = 1;

% Define the function to add artificial delay δ to Yr
Yr_delayed      = @(delta) -abs(e_k.' * (F * (Yr_tilde .* exp(-1j * 2 * pi * delta * (0:Nsc-1)' * Df)))) / y0_hat;

% Perform the optimization over δ in the range [0, Delta_t]
delta_hat       = fminunc(Yr_delayed, 0, options);

% Estimate tau_r
tau_r_est_coarse = (k_hat - 1) * Delta_t;
tau_r_est_fine   = tau_r_est_coarse - delta_hat;




% % OFDM symbol generation
% data    = 2*randi([0 1],1, Nsym*Nfft)-1; % BPSK data
% Tx      = zeros(1,Nsym*(Nfft+Ncp));
% OFDMsym = zeros(1,Nfft);  
% 
% for sym = 1:Nsym
%     OFDMsym                                 = ifft(data(Nfft*(sym-1)+1:(Nfft*sym)),Nfft)*sqrt(Nfft);
%     Tx((Nfft+Ncp)*(sym-1)+1:(Nfft+Ncp)*sym) = [OFDMsym(Nfft-Ncp+1:Nfft) OFDMsym];
% end
% 
% %AWGN channel
% noise = sqrt(snr/2)*(randn(1,Nsym*(Nfft+Ncp))+1i*randn(1,Nsym*(Nfft+Ncp)));
% Rx    = exp(1i*2*pi*FreqOffset*(0:length(Tx)-1)./Nfft).*Tx + noise; 


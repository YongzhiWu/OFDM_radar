clc;
clear;

c0 = physconst('LightSpeed');
fc = 30e9; % carrier frequency
lambda = c0 / fc; % wavelength
N = 64; % number of subcarriers
M = 16; % number of symbols
delta_f = 15e3 * 2^6; % subcarrier spacing
T = 1 / delta_f; % symbol duration
Tcp = T / 4; % cyclic prefix duration
Ts = T + Tcp; % total symbol duration

qam = 4; % 4-QAM modulation

% Transmit data
data = randi([0 qam - 1], N, M);
TxData = qammod(data, qam, 'UnitAveragePower', true);

% Target parameters
target_pos = 50; % target distance
target_delay = range2time(target_pos, c0);
target_speed = 20; % target velocity
target_dop = speed2dop(2 * target_speed, lambda);

SNR_dB = 5;
SNR = 10.^(SNR_dB/10);

% Received frequency-domain signal over the radar channel
RxData = zeros(size(TxData));
for kSubcarrier = 1:N
    for mSymbol = 1:M
        RxData(kSubcarrier, mSymbol) = sum(sqrt(SNR) * TxData(kSubcarrier, mSymbol) .* exp(-1j * 2 * pi * fc * target_delay)  .* exp(1j * 2 * pi * mSymbol *...
            Ts .* target_dop) .* exp(-1j * 2 * pi * kSubcarrier .* target_delay *...
            delta_f) ) + sqrt(1/2)* (randn() +1j * randn());
    end
end

% Remove the information of transmit data
dividerArray = RxData ./ TxData;

%% MUSIC algorithm
nTargets = 1;
Rxxd = dividerArray * dividerArray' / M;
[Vd, distanceEigen] = eig(Rxxd);
distanceEigenDiag = diag(distanceEigen)';
[distanceEigenDiag, eigenMark] = sort(distanceEigenDiag);
distanceEigenDiag = fliplr(distanceEigenDiag);
Vd = fliplr(Vd(:, eigenMark));
omegaDistance = 0 : pi / 100 : 2 * pi;
distanceEigenMatNoise = Vd(:, nTargets + 1: end);
SP = zeros(1, length(omegaDistance));
nIndex = 0:1:N - 1;
for index = 1:length(omegaDistance)
    omegaVector = exp(-1i * nIndex * omegaDistance(index)).';
    SP(index) = (omegaVector' * omegaVector ) / (omegaVector' * (distanceEigenMatNoise * distanceEigenMatNoise') * omegaVector);
end
SP = abs(SP);
SPmax = max(SP);
SP_dB = 10*log10(SP/SPmax);
distanceIndex = omegaDistance * c0 / (2 * pi * 2 * delta_f);
plot(distanceIndex, SP_dB);

%% Periodogram/FFT-based estimation
NPer = 16 * N;
normalizedPower = abs(ifft(dividerArray, NPer, 1));
normalizedPower_dB = 10 * log10(normalizedPower);
mean_normalizedPower = mean(normalizedPower, 2);
mean_normalizedPower = mean_normalizedPower / max(mean_normalizedPower);
mean_normalizedPower_dB = 10 * log10(mean_normalizedPower);
[~, rangeEstimation] = max(mean_normalizedPower_dB);
rangeIndex = 0:1:(NPer-1);
rangeIndex = rangeIndex * c0 / (2 * delta_f * NPer);
distanceE = rangeEstimation * c0 / (2 * delta_f * NPer);
hold on;
plot(rangeIndex, mean_normalizedPower_dB);
grid on;
xlabel('Range [m]');
ylabel('Normalized Range Profile [dB]');
legend('MUSIC', 'periodogram')

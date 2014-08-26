%This Matlab script can be used to generate Figure 5, in the article:
%
%Emil Björnson, Jakob Hoydis, Marios Kountouris, Mérouane Debbah, “Massive
%MIMO Systems with Non-Ideal Hardware: Energy Efficiency, Estimation, and
%Capacity Limits,” To appear in IEEE Transactions on Information Theory.
%
%Download article: http://arxiv.org/pdf/1307.2584
%
%This is version 1.0 (Last edited: 2014-08-26)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


%Initialization
close all;
clear all;


%%Simulation parameters

%Range of number of BS antennas
N = [1:6 8:2:100 105:3:200 200];
Nmax = max(N); %Maximal number of antennas in simulation

%Number of different channel covariance models
nbrOfCovarianceModels = 4;

%Range of SNRs in simulation (we have normalized sigma2 to 1)
SNRdB = [5 30]; %In decibel scale
SNR = 10.^(SNRdB/10); %In linear scale

%Define the level of hardware impairments at the transmitter and receiver
kappatUE = 0.05^2;
kapparBS = 0.05^2;


%Placeholders for storing simulation results (normalized MSE for LMMSE
%estimator in Theorem 1)
normalizedMSE_LMMSE = zeros(length(N),nbrOfCovarianceModels,length(SNR));


%Go through different number of antennas
for n = 1:length(N)
    
    %Output the progress of the simulation
    disp(['Antennas: ' num2str(n) '/' num2str(length(N))]);
    
    %Go through SNR values
    for m = 1:length(SNR)
        
        %Compute the pilot signal for the given SNR value
        d = sqrt(SNR(m));
        
        %Go through the different channel covariance models
        for c = 1:nbrOfCovarianceModels
            
            %Generate channel covariance matrix
            if c == 1 %One-ring model with 10 degrees angular spread
                R = functionOneRingModel(N(n),10);
            elseif c == 2 %One-ring model with 20 degrees angular spread
                R = functionOneRingModel(N(n),20);
            elseif c == 3 %Uncorrelated covariance matrix
                R = eye(N(n));
            elseif c == 4 %Exponential correlation model
                correlationFactor = 0.7;
                R = toeplitz(correlationFactor.^(0:N(n)-1));
            end
            
            %Compute matrix A in the LMMSE estimator (see Eq. (9))
            A_LMMSE = conj(d) * R(1:N(n),1:N(n)) / (abs(d)^2*(1+kappatUE)*R(1:N(n),1:N(n)) + abs(d)^2*kapparBS*diag(diag(R(1:N(n),1:N(n))))+eye(N(n)));
            
            %Compute the MSE as in Eq. (12) and normalize by tr(R)
            normalizedMSE_LMMSE(n,c,m) = real(1 - trace(d*  A_LMMSE *R(1:N(n),1:N(n)))/trace(R(1:N(n),1:N(n))));
            
        end
        
    end
    
end



%Plot Figure 5 from the paper
figure; hold on; box on;

for m = 1:length(SNR)
    plot(N,normalizedMSE_LMMSE(:,3,m),'k','LineWidth',1);
    plot(N,normalizedMSE_LMMSE(:,4,m),'r:','LineWidth',1);
    plot(N,normalizedMSE_LMMSE(:,2,m),'b-.','LineWidth',1);
    plot(N,normalizedMSE_LMMSE(:,1,m),'k--','LineWidth',1);
end

set(gca,'YScale','log');
xlabel('Number of Base Station Antennas (N)');
ylabel('Relative Estimation Error per Antenna');
legend('Case 1: Uncorrelated','Case 2: Exponential Mod. r=0.7','Case 3: One-Ring, 20 degrees','Case 4: One-Ring, 10 degrees','Location','SouthWest');

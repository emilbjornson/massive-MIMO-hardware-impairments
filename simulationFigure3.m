%This Matlab script can be used to generate Figure 3, in the article:
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

N = 50; %Number of BS antennas

%Compute normalized channel covariance matrix R according to the exponential
%correlaton model in Eq. (17).
correlationFactor = 0.7;
R = toeplitz(correlationFactor.^(0:N-1));

%Define the range of level of hardware impairments values (same at
%transmitter and receiver for simplicity)
kappaValues = [0 0.05^2 0.1^2 0.15^2];

%Range of SNRs in simulation (we have normalized sigma2 to 1)
SNRdB = 0:5:40; %In decibel scale
SNR = 10.^(SNRdB/10); %In linear scale


%Placeholders for storing simulation results
normalizedMSE_LMMSE = zeros(length(SNR),length(kappaValues)); %Normalized MSE for LMMSE estimator in Theorem 1
normalizedMSE_impairment_ignore = zeros(length(SNR),length(kappaValues)); %Normalized MSE for Conventional estimator for kappa=0


%Go through SNR values
for m = 1:length(SNR)
    
    %Compute the pilot signal for the given SNR value
    d = sqrt(SNR(m));
    
    %Go through all level of impairments values
    for impairIndex = 1:length(kappaValues)
        
        %Extract current level of impairments at transmitter and receiver
        kappatUE = kappaValues(impairIndex);
        kapparBS = kappaValues(impairIndex);
        
        
        %Compute matrix A in the LMMSE estimator (see Eq. (9))
        A_LMMSE = conj(d) * R(1:N,1:N) / (abs(d)^2*(1+kappatUE)*R(1:N,1:N) + abs(d)^2*kapparBS*diag(diag(R(1:N,1:N)))+eye(N));
        
        %Compute the MSE as in Eq. (12) and normalize by tr(R)
        normalizedMSE_LMMSE(m,impairIndex) = 1 - trace(d*  A_LMMSE *R(1:N,1:N))/trace(R(1:N,1:N));
        
        
        %Compute matrix A for impairment-ignoring case by setting kappat=kappar=0 in Eq. (9)
        A_ignore = conj(d) * R(1:N,1:N) / (abs(d)^2*R(1:N,1:N) +eye(N));
        
        %Compute the MSE as in Eq. (12), by setting kappat=kappar=0 , and normalize by tr(R) 
        normalizedMSE_impairment_ignore(m,impairIndex) = trace(R(1:N,1:N)- d*A_ignore*R(1:N,1:N) - d'*R(1:N,1:N)*A_ignore' + A_ignore* (abs(d)^2*(1+kappatUE)*R(1:N,1:N) + abs(d)^2*kapparBS*diag(diag(R(1:N,1:N)))+eye(N)) *A_ignore')/trace(R(1:N,1:N));
        
    end
    
end



%Compute the asymptotic limits (similar to Corollary 3, but for general R)
asymptoticLimits = zeros(length(kappaValues),1);

for impairIndex = 1:length(kappaValues) %Go through all level of impairments
    
    %Extract current level of impairments at transmitter and receiver
    kappatUE = kappaValues(impairIndex);
    kapparBS = kappaValues(impairIndex);
    
    %Compute matrix A in the LMMSE estimator (see Eq. (9)) for large SNR
    %when noise contribution can be ignored.
    A_LMMSE = conj(d) * R(1:N,1:N) / (abs(d)^2*(1+kappatUE)*R(1:N,1:N) + abs(d)^2*kapparBS*diag(diag(R(1:N,1:N))));
    
    
    asymptoticLimits(impairIndex) = 1 - trace(d*  A_LMMSE *R(1:N,1:N))/trace(R(1:N,1:N));
    
    
end



%Plot Figure 3 from the paper
figure; hold on; box on;

for impairIndex = 1:length(kappaValues)
    
    plot(SNRdB,asymptoticLimits(impairIndex)*ones(size(SNRdB)),'r:','LineWidth',1);
    
    plot(SNRdB,normalizedMSE_LMMSE(:,impairIndex),'ko-','LineWidth',1);
    
    plot(SNRdB,normalizedMSE_impairment_ignore(:,impairIndex),'bd-.','LineWidth',1);
    
end

set(gca,'YScale','log'); %Logarithmic vertical scaling
xlabel('Average SNR [dB]');
ylabel('Relative Estimation Error per Antenna');
legend('Conventional Impairment-Ignoring','LMMSE Estimator','Error Floors','Location','SouthWest');

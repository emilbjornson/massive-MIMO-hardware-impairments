function R = functionOneRingModel(N,angularSpread)
%This is an implementation of the channel covariance matrix with the
%one-ring model. The implementation is based on Eq. (57) in the paper:
%
%A. Adhikary, J. Nam, J.-Y. Ahn, and G. Caire, “Joint spatial division and
%multiplexing—the large-scale array regime,” IEEE Trans. Inf. Theory,
%vol. 59, no. 10, pp. 6441–6463, 2013.
%
%This is used in the article:
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
%
%INPUT:
%N             = Number of antennas
%angularSpread = Angular spread around the main angle of arrival
%
%OUTPUT:
%R             = N x N channel covariance matrix


%Approximated angular spread
Delta = angularSpread*pi/180; 

%Half a wavelength distance
D = 1/2; 

%Angle of arrival (30 degrees)
theta = pi/6; 

%The covariance matrix has the Toeplitz structure, so we only need to
%compute the first row.
firstRow = zeros(N,1);

%Go through all columns in the first row
for col = 1:N
    
    %Distance from the first antenna
    distance = col-1;
    
    %Define integrand of Eq. (57) in [42]
    F = @(alpha)exp(-1i*2*pi*D*distance*sin(alpha+theta))/(2*Delta);
    
    %Compute the integral as in [42]
    firstRow(col) = integral(F,-Delta,Delta);
    
end

%Compute the covarince matrix by utilizing the Toeplitz structure
R = toeplitz(firstRow);

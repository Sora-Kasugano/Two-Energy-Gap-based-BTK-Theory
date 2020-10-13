function  [S, DI1] = mBTK2(aa,V)
% The matlab program is coded by github user Sora-Kasugano.
% All rights reserved.
%
% This program is intended for the experimental fit of modified BTK theory
% considering the Two-Energy Gap assumption. Related work stemmed from the
% paper: PhysRevB.63.104510. The following codes are based from this paper.
% 
% Tips of this program:
% 1: ''V'' stands for a group of DC bias voltage data. By default, the unit of V is mV.
% 2: ''aa'' stands for a row vector, which contains the following fit parameters.
% a= aa(1); b= aa(2); Delta1= aa(3); Delta2= aa(4); Gama=aa(5); Z= aa(6); P= aa(7); npanel= aa(8); T= aa(9);
%
%      a ,b          = lower and upper limits of the integral.
%      Delta         = Superconducting energy gap. In Two-Energy Gap based
%                      modified BTK theory, Delta1 is no more than
%                      Delta2 and their units are mV.
%      Gama          = Inelastic scattering factor with its unit mV.
%      Z             = Barrier factor
%      P             = Spin Polarization
%      T             = Experimental Temperature (K)
%      npanel        = number of panels to use in the integration, with
%                      Total number of nodes = 2*panel + 1
% 
a= aa(1); b= aa(2); Delta1= aa(3); Delta2= aa(4); Gama=aa(5); Z= aa(6); P= aa(7); npanel= aa(8); T= aa(9); 
% Every V gets a value of Andreev Reflection Current. 
for m=1:length(V)
    n = 2 * npanel + 1;                                          % total number of nodes
    h = (b-a)/(n-1);   
    % stepsize
    E = a:h:b;                                                         % divide the interval
    u01 = 0.5 * (1 + real(sqrt(((E - i * Gama).^2 - Delta1^2)./(E - i * Gama).^2)));             % coefficient of BCS u01^2
    v01 = 0.5 * (1 - real(sqrt(((E - i * Gama).^2 - Delta1^2)./(E - i * Gama).^2)));              % coefficient of BCS v01^2
    u02 = 0.5 * (1 + real(sqrt(((E - i * Gama).^2 - Delta2^2)./(E - i * Gama).^2)));             % coefficient of BCS u02^2
    v02 = 0.5 * (1 - real(sqrt(((E - i * Gama).^2 - Delta2^2)./(E - i * Gama).^2)));              % coefficient of BCS v02^2
    gama1 = (u01 + Z^2 * (u01 - v01)).^2;
    gama2 = u01.* v01 + (u02 - v02).* (u02 + Z^2 + (u02 - u01) * Z^2 * (1+Z^2));
    gama3 = (u02 - v02).* (u02 + Z^2 + (u02 - u01) * Z^2 * (1+Z^2));                           % coefficient of BTK
    % define the coefficient of BTK theory A, B
    j = 1;
    while j <= length(E)
        if abs(E(j)-i*Gama) <= Delta1
            AN(j) = real(Delta1^2/((E(j) - i * Gama)^2 + (Delta1^2 - (E(j) - i * Gama)^2) * (1 + 2 * Z^2)^2));
            BN(j) = 1 - AN(j);
            BP(j) = 1;
        elseif (abs(E(j)-i*Gama) )<= Delta2
            AN(j) = u01(j) * v01(j)/gama1(j);
            BN(j) = 1 - AN(j);
            BP(j) = 1;
        else
            AN(j) = u01(j) * v01(j)/gama2(j);
            BN(j) = (u02(j) - v02(j))^2 * Z^2 * (1 + Z^2)/gama2(j);
            BP(j) = (u02(j) - v02(j))^2 * Z^2 * (1 + Z^2)/gama3(j);
        end
        j = j + 1;
    end
    f = 1./(1+exp(11.5942*(E-V(m))/T)) - 1./(1+exp(11.5942*E/T));                              
    z = ((1-P) * (1 + AN - BN) + P * (1 - BP)).* f;
    I(m) = (h/3)*(z(1)+4*sum(z(2:2:n-1))+2*sum(z(3:2:n-2))+z(n));
end
% Claculate the normalized differential conductance and output them with
% the form [S DI1], where S represents the DC bias voltage, and DI1 the
% normalized differential conductance. 
S=V;
S(1)=[]; % Ensure the same dimention of normalized differential conductance and DC bias voltage.
DI=diff(I);
DI1=DI/DI(length(V)-1); 
plot(S,DI1)
% Since S and DI1 are in a same dimension, it's easy to plot the result fit of
% Two-Energy Gap based modified BTK theory.
end

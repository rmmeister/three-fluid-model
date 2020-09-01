function [ F, GammaDL, GammaLD, dropletMeanDiameter, RGD, RLW, RGL, RDL ] = calcResiduals(  Y, alphaGvG )
%CALCRESIDUALS Calulate the residuals vector, F.
%   Residual terms are obtained using well data and vector Y containing
%   parameters to be defined.
gravityAcceleration = 9.8; % m/s2

% INPUT DATA
pipeDiameter = .032; % m
kLin = 80; % Liquid mass flow-rate, Kg/s/m2
kGin = alphaGvG; % Gas mass flow-rate, Kg/s/m2
deviationAngle = degtorad(90); % deviation angle

% pvt
gasViscosity = 1.78e-5; % Kg/m/s
liquidViscosity = 9.98e-4; % Kg/m/s
interfacialTension = .073; % N/m
gasDensity = 1.725; % Kg/m3
liquidDensity = 998; % Kg/m3

% Parameters to be found
alphaG = Y(1); alphaL = Y(2); alphaD = Y(3); vG = Y(4); vL = Y(5); vD = Y(6);
gasCoreDiameter = pipeDiameter*(1 - alphaL)^(1/2);

% calculate frictional CLOSURES
[ RGD, RLW, RGL, RDL, GammaDL, GammaLD, dropletMeanDiameter  ] = calcTFM( alphaG, alphaL, alphaD, vG, vL, vD, gasDensity,...
    liquidDensity, gasViscosity, liquidViscosity, pipeDiameter, gasCoreDiameter, deviationAngle, interfacialTension );

% relevant parameters for residuals
vL_D = vL ;
vD_L = vD;
FG = -RGL - RGD - alphaG*gasDensity*gravityAcceleration*sin(deviationAngle);
FL = -vL_D*GammaLD + vD_L*GammaDL + RGL + RDL - RLW - alphaL*liquidDensity*gravityAcceleration*sin(deviationAngle);
FD = vL_D*GammaLD - vD_L*GammaDL - RDL + RGD - alphaD*liquidDensity*gravityAcceleration*sin(deviationAngle);

% obtain RESIDUAL terms
R1 = alphaG*vG - kGin;
R2 = (alphaL*vL + alphaD*vD)*liquidDensity - kLin;
R3 = alphaG + alphaL + alphaD - 1;
R4 = GammaLD - GammaDL;
R5 = FG/alphaG - FL/alphaL;
R6 = FG/alphaG - FD/alphaD;

% vector of RESIDUALS
F = [R1; R2; R3; R4; R5; R6];
end


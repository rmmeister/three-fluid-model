function [ RGD, RLW, RGL, RDLi, GammaDL, GammaLD, dD  ] = calcTFM( alphaG, alphaL, alphaD, vG,...
    vL, vD, gasDensity, liquidDensity, gasViscosity, liquidViscosity, d, di, deviationAngle, interfacialTension )
%CALCTPM Calculate three-fluid model properties
%

ks = 0; % surface roughness
gravityAcceleration = 9.8;

% determine friction factor
deltaL = (d - di)/2;
ReGL = (vG-vL)*di*gasDensity/gasViscosity;
% Haaland's explicit approximation of Colebrook and White's equation
f = ( -1.8*log10(6.9/ReGL + (max([0 ks-deltaL])/3.7/d)^1.11) )^-2;
fGL = f + 6*deltaL/d;
ReLW = alphaL*abs(vL)*d*liquidDensity/liquidViscosity;
taoGL = alphaG/(alphaG + alphaD)*fGL*gasDensity/8*(vG-vL)*abs(vG-vL);

% Calculate mean droplet diameter
dDdvs = getMaxDdVelocity ( 1e-4, liquidDensity, gasDensity,interfacialTension, liquidViscosity, vG, vL);
dDturbs = getMaxDdTurb ( 1e-4, liquidDensity, gasDensity, interfacialTension, gasViscosity, vG, vD, taoGL, di);
dDs = min([dDdvs dDturbs]);
dD = .25*dDs;

% particle relaxation time (droplet's dynamic response time)
ReDG = abs(vG - vD)*dD*gasDensity/gasViscosity;
if ReDG < .2
    CD = 24/ReDG;
else
    CD = 24/ReDG*(1+.27*ReDG)^.43 + .47*( 1 - exp(-0.04*ReDG^.38) );
end
CDm = .5;
tDR = 4*( liquidDensity + CDm*gasDensity )*dD/( 3*gasDensity*CD*abs(vG-vD) );

% eddy droplet interaction time model
tE = .1*di/ ( (abs(taoGL)/gasDensity)^(1/2) );
tLag = .4*tE;
St = tDR/tE;
vGDstar = abs(vG-vD)/ ( abs(taoGL)/gasDensity )^(1/2);
Ca = ( 1 + vGDstar^2 + 2*vGDstar/sqrt(3) )^(1/2);
tDi0 = tLag * (12*Ca + 6*Ca^2 + 2 ) / ( 5*Ca*(1+Ca)^2 );
tDiInf = tLag* 6*( 2+vGDstar ) / ( 5* (1+vGDstar)^2 );
fSt = St/(1+St) - 5*St^2/ ( 4*(1+St)^2*(2+St) );
tDi = tDi0 + (tDiInf - tDi0)*fSt;

% Droplet-liquid film friction model (continuum droplet)
C0 = ( 1 + CDm )*(gasDensity/liquidDensity) / ( 1 + CDm*gasDensity/liquidDensity );
CGR = (1 + C0*tDR/tDi ) / ( 1 + tDR/tDi );
taoDL = alphaD*liquidDensity*(vD-vL)^2*CGR*taoGL / ( alphaG*gasDensity*(vG-vL)^2 );

% Haaland's explicit approximation of Colebrook and White's equation
f = ( -1.8*log10(6.9/ReLW + (ks/3.7/d)^1.11) )^-2;
% friction between liquid film and wall
G = (liquidDensity - alphaG*gasDensity - alphaD*liquidDensity)*deltaL/(taoGL+taoDL)*gravityAcceleration*abs(sin(deviationAngle));
kf = (1+G)/(1+2/3*G);
kd = (1 + exp( (ReLW - 1600)/100 ) )^(-1);
fLW = (kd*kf + (1-kd))*f;

% friction closures
RGD = 3*alphaD*gasDensity*CD/4/dD*(vG-vD)*abs(vG-vD);
RLW = fLW*liquidDensity/2/d*vL*abs(vL);
RGL = di/d^2*alphaG/(alphaG+alphaD)*fGL*gasDensity/2*(vG-vL)*abs(vG-vL);
RDLi = 4*di/d^2*taoDL;

% Droplet deposition
Cr = 5/81;
vGr2 = ( tDi*Cr*di / (tLag*2*tDR) * (abs(taoGL)/gasDensity)^(1/2) ) / ...
    ( 1 + ( tDi*16*Cr*di*tDR / (tLag*pi) * (abs(taoGL)/gasDensity)^(1/2) * (alphaD/5/dD + 1/di)^2 )^(1/3) );
WeLGD = ( taoGL + taoDL ) * deltaL / interfacialTension;
X = 1 - exp ( -3*WeLGD );
gammaNum = 32*alphaD*liquidDensity*tDR*vGr2/d^2;
gammaDenom = 1 + ( (1-X)/(1+X)*sqrt(2/pi) + 4*(1+X)/(1-X)*sqrt(pi/2) ) * ( 2*tDR*sqrt(vGr2) ) / di;
GammaDL = gammaNum / gammaDenom;

% Liquid film entrainment
ReLWs = 160; % or 160
WeLGDs = ( 7e-6 + 4e-4*( ReLW - ReLWs )^(-.8) )*ReLW;
if WeLGD > WeLGDs
    GammaLD = .023*4*di/d^2*sqrt( liquidDensity* (taoGL + taoDL) )* ( WeLGD - WeLGDs );
elseif WeLGD <= WeLGDs
    GammaLD = 0;
end

end

function dDdvs = getMaxDdVelocity ( dDdvs, rhoL, rhoG,sigmaLG, muL, vG, vL)
maxIter = 50;
velError = inf;
iter = 0;
while velError > 1e-12
    dD = dDdvs;
    WeDGs = 12 + 18* ( rhoL*sigmaLG*dD/muL^2 )^(-.37);
    dvGD = vG - vL;
    % WeDG = dD*rhoG*dvGD^2/sigmaLG;
    dDdvsOld = dDdvs;
    dDdvs = WeDGs*sigmaLG/rhoG/dvGD^2;
    velError = abs(dDdvs - dDdvsOld);
    iter = iter + 1;
    if iter > maxIter
        fprintf('\nMaximum Iteration Reached on Droplet Diameter in Velocity Mode\n');
        fprintf('the Error is: %g\n', velError);
        break
    end
end
end

function dDturbs = getMaxDdTurb ( dDturbs, rhoL, rhoG, sigmaLG, muG, vG, vD, taoGL, di)
maxIter = 100;
iter = 0;
turbsError = inf;
while turbsError > 1e-12
    dD = dDturbs;
    n = 20;
    % particle relaxation time (droplet's dynamic response time)
    ReDG = abs(vG - vD)*dD*rhoG/muG;
    if ReDG < .2
        CD = 24/ReDG;
    else
        CD = 24/ReDG*(1+.27*ReDG)^.43 + .47*( 1 - exp(-0.04*ReDG^.38) );
    end
    CDm = .5;
    tDR = 4*( rhoL + CDm*rhoG )*dD/( 3*rhoG*CD*abs(vG-vD) );
    
    Cmu = .09;
    kG = (taoGL/rhoG)/sqrt(Cmu);
    epsilonG = 25*((taoGL/rhoG)^(3/2))/di;
    
    % eddy droplet interaction time model
    tE = .1*di/ ( (abs(taoGL)/rhoG)^(1/2) );
    tLag = .4*tE;
    St = tDR/tE;
    vGDstar = abs(vG-vD)/ ( abs(taoGL)/rhoG )^(1/2);
    Ca = ( 1 + vGDstar^2 + 2*vGDstar/sqrt(3) )^(1/2);
    tDi0 = tLag * (12*Ca + 6*Ca^2 + 2 ) / ( 5*Ca*(1+Ca)^2 );
    tDiInf = tLag* 6*( 2+vGDstar ) / ( 5* (1+vGDstar)^2 );
    fSt = St/(1+St) - 5*St^2/ ( 4*((1+St)^2)*(2+St) );
    tDi = tDi0 + (tDiInf - tDi0)*fSt;
    
    C0 = ( 1 + CDm )*rhoG/rhoL / ( 1 + CDm*rhoG/rhoL );
    SG = ( (15*muG/(epsilonG*dD^2*rhoG))^n + (1/( 2*epsilonG^(2/3)*dD^(2/3) ))^n + (3/4/kG)^n )^(-1/n);
    tdDi = ( 1 / ( (5*muG/epsilonG/rhoG)^(n/2) + (.3*dD^(2/3)/epsilonG^(1/3))^n ) + (1/tDi)^n )^(-1/n);
    
    CGR2 = (1 + C0*tDR/tdDi ) / ( 1 + tDR/tdDi );
    
    WeDGturbs = 3 / ( CGR2*(1 + 2*rhoL/rhoG) );
    dDturbsOld = dDturbs;
    dDturbs = WeDGturbs * sigmaLG/rhoG/SG;
    turbsError = abs(dDturbsOld - dDturbs);
    iter = iter + 1;
    if iter > maxIter
        fprintf('\nMaximum Iteration Reached on Droplet Diameter in Turbulent Mode\n');
        fprintf('the Error is: %g\n', turbsError);
        break
    end
end
end
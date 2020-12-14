function [r0, v0, OE0, rf, vf, OEf] = OrbitComp(phi, lambda, rho, beta...
    ,sigma,rho_dot, beta_dot, sigma_dot, TOF)
     
    %{
    %testers
    phi = 32.248814; %latitude of ground station [deg]
    lambda = -74.99; %longitude of ground station [deg]
    rho = 822; %range [km]
    beta = 18.0912; %azimuth angle [deg]
    sigma = 61.7066; %elevation [deg]
    rho_dot = 3.48499169; %range rate [km/s]
    beta_dot = 0.269604966; %azimuth rate [deg/s]
    sigma_dot = -0.4321605433; %elevation rate [deg/s]
    TOF = 32*60; %10*60*60; %time of flight [s]
    %}
    %rotation of earth and mu
    omegaE = [0; 0; 7.2921*10^-5]; %[rad/s]
    mu = 398600.44; 
    RE = 6378.1366; %[km]
    
    %convert degrees to radians 
    phi = phi*pi/180;
    lambda = lambda*pi/180;
    beta = beta*pi/180;
    sigma = sigma*pi/180;
    beta_dot = beta_dot*pi/180;
    sigma_dot = sigma_dot*pi/180;
    
    %find the SEZ components of rho and rho_dot
    [rhoSEZ,rhodotSEZ] = SEZrhofun(beta,sigma,rho,beta_dot,sigma_dot,rho_dot);
    
    %find coordinates in the ECI frame 
    [rhoECI] = SEZECISEZ(lambda,phi,rhoSEZ,'ECI');
    [rhodotECI] = SEZECISEZ(lambda,phi,rhodotSEZ,'ECI');
    
    %find the r and v vectors in ECI frame 
    %rsite assumed to be at sea level
    rsite = [0; 0; RE];
    [rsiteECI] = SEZECISEZ(lambda,phi,rsite,'ECI');
    [r0,v0] = rhotoRV(rsiteECI,rhoECI,rhodotECI,omegaE);
    
    %find orbital elements from r and v
    [a,e,i,BOmegai,lomegai,fi,h,e_hat,n_hat,h_hat,ti]  = RVtoOE(r0,v0,mu);
    
    %initial orbital elements vector
    OE0 = [a; e; i*180/pi; BOmegai*180/pi; lomegai*180/pi; fi*180/pi];
    
    %add time to find the f start with new M
    n = sqrt(mu/a^3);
    M = n*(ti+TOF);
    
    %newtons method to find the new f
    [ff,E,F] = NewtonMethod(M,0.01,e);
    
    %adjust OE for pertibations 
    J2 = 1.08264*10^-3;
    BOmega_dot = -((3*n*J2)/(2*(1-e^2)^2)*(RE/a)^2*cos(i));
    lomega_dot = -((3*n*J2)/(4*(1-e^2)^2)*(RE/a)^2*(5*(cos(i))^2-1));
    
    BOmegaf = BOmegai+(BOmega_dot*TOF);
    lomegaf = lomegai+(lomega_dot*TOF);
    
    %final orbital elements
    OEf = [a; e; i*180/pi; BOmegaf*180/pi; lomegaf*180/pi; ff*180/pi];
    
    %final r and v
    [rf,vf]= OEtoRV(a,e,i,BOmegaf,lomegaf,ff,mu);
    
end


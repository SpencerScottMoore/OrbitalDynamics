function [rhoSEZ,rhodotSEZ] = SEZrhofun(beta,sigma,rho,betadot,sigmadot,rhodot)
    %finds components of rho in SEZ
    rhos = -rho*cos(sigma)*cos(beta);
    rhoe = rho*cos(sigma)*sin(beta);
    rhoz = rho*sin(sigma);
    
    %combines them all
    rhoSEZ = [rhos; rhoe; rhoz];
    
    %finds components of rhodot in SEZ frame
    rhodots = -rhodot*cos(sigma)*cos(beta)+rho*sigmadot*sin(sigma)...
        *cos(beta)+rho*betadot*cos(sigma)*sin(beta);
    rhodote = rhodot*cos(sigma)*sin(beta)-rho*sigmadot*sin(sigma)...
        *sin(beta)+rho*betadot*cos(sigma)*cos(beta);
    rhodotz = rhodot*sin(sigma)+rho*sigmadot*cos(sigma);
    
    %combines them
    rhodotSEZ = [rhodots; rhodote; rhodotz];

end


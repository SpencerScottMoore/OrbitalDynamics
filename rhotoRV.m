function [r_vec,v_vec] = rhotoRV(rsite,rho,rho_dot,omegaE)
    %coverts ground track vectors to state vectors
    
    %makes the conversions
    r_vec = rsite+rho;
    v_vec = cross(omegaE,rsite)+rho_dot+cross(omegaE,rho);
    
end


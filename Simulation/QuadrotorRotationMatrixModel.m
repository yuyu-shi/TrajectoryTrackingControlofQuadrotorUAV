function [pNext, vNext, RNext, omegaNext] = QuadrotorRotationMatrixModel(p, v, R, omega, f, tau, df, dtau)

    global dt m J JInv g e3;    

    pNext = dt * v + p;
    vNext = dt * (g * e3 - f * R(:, 3)/ m  + df / m) + v;
    RNext = R * expm(R3so3(dt * omega));
    omegaNext = dt * JInv * ( - cross(omega, J * omega) + tau + dtau) + omega;
    
end
function [pNext, vNext, ThetaNext, omegaNext] = QuadrotorEulerAngleModel(p, v, Theta, omega, f, tau, df, dtau)

    global dt m J JInv g e3;
    phi = Theta(1);
    theta = Theta(2);
    psi = Theta(3);
    
    W = [1, tan(theta) * sin(phi), tan(theta) * cos(phi);
         0, cos(phi), -sin(phi);
         0, sin(phi) * sec(theta), cos(phi) * sec(theta)];
     
    R = [cos(theta) * cos(psi), cos(psi) * sin(theta) * sin(phi) - sin(psi) * cos(phi), cos(psi) * sin(theta) * cos(phi) + sin(psi) * sin(phi);
         cos(theta) * sin(psi), sin(psi) * sin(theta) * sin(phi) + cos(psi) * cos(phi), sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi);
         -sin(theta), sin(phi) * cos(theta), cos(phi) * cos(theta)];
    pNext = dt * v + p;
    vNext = dt * (g * e3 - f / m * R * e3 + 1 / m * df) + v;
    ThetaNext = dt * W * omega + Theta;
    omegaNext = dt * JInv * (-cross(omega, J * omega) + tau + dtau) + omega;
end
function [pNext, vNext, QNext, omegaNext] = QuadrotorQuaternionModel(p, v, Q, omega, f, tau, df, dtau)

    global dt m J JInv g e3;
    
    R = eye(3) + 2 * R3so3(Q(2:4))^2 + 2 * Q(1) * R3so3(Q(2:4));
    
    pNext = dt * v + p;
    vNext = dt * (g * e3 - f / m * R * e3 + 1 / m * df) + v;
    if norm(omega) == 0
        QNext = Q;
    else
        omegatheta = norm(omega);
        omegaa = omega / omegatheta;
        DeltaQ = [cos(omegatheta/2*dt);
                  sin(omegatheta/2*dt) * omegaa];
        DeltaQ = DeltaQ / norm(DeltaQ);
        QNext = Quatprod(Q, DeltaQ, 0);
    end
    omegaNext = dt * JInv * ( - cross(omega, J * omega) + tau + dtau) + omega;    
end

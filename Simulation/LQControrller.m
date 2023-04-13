function [f, tau, state] = LQControrller(t, p, v, pose, omega, pd, psid, attFlag, enumR, enumA)
global g e3 m J;
    persistent Kp Kg;
    if isempty(Kp)
        I = eye(3);
        O = zeros(3,3);
        Atilde = [O I;
                  O O];
        Btilde = [O;
                  I];
        Qp = 10 * eye(6);
        Rp = 0.01 * eye(3);
        Qg = 5 * eye(6);
        Rg = 0.01 * eye(3);
        Kp = lqr(Atilde, Btilde, Qp, Rp);
        Kg = lqr(Atilde, Btilde, Qg, Rg);
    end
    ptilde = p - pd(:,1);
    vtilde = v - pd(:,2);
    xptilde = [ptilde;
               vtilde];

    FdLump = -Kp * xptilde;
    Fd = -m * (FdLump - g * e3 + pd(:,3)); 
%     f = Fd' * (ToSO3(pose, attFlag) * e3);
    
    [gtilde, G, GInv, omegatilde, A, omegadbar, data] = GibbsVector(t, pose, Fd, psid, omega, enumR, enumA, attFlag);
    Fd = data{1};
    f = Fd' * (ToSO3(pose, attFlag) * e3);
    
    gtilde1d = G * omegatilde;
    xgtilde = [gtilde;
               gtilde1d];
    tauLump = -Kg * xgtilde;
    tau = J * A' * (GInv * (tauLump - 2 * gtilde' * gtilde1d * gtilde1d / (1 + gtilde' * gtilde)) - omegadbar) + cross(omega, J * omega);
    
    varthetatilde = 2 * atan(norm(gtilde));
    state{1} = ptilde;
    state{2} = varthetatilde;
end
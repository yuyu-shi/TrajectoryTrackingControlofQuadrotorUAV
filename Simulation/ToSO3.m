function R = ToSO3(pose, Flag)
    if Flag == 'vartheta'
        R = expm(R3so3(pose));
    elseif Flag == 'Q'
        R = eye(3) + 2 * R3so3(pose(2:4))^2 + 2 * pose(1) * R3so3(pose(2:4));
    elseif Flag == 'Theta'
        phi = pose(1);
        theta = pose(2);
        psi = pose(3);
        R = [cos(theta) * cos(psi), cos(psi) * sin(theta) * sin(phi) - sin(psi) * cos(phi), cos(psi) * sin(theta) * cos(phi) + sin(psi) * sin(phi);
         cos(theta) * sin(psi), sin(psi) * sin(theta) * sin(phi) + cos(psi) * cos(phi), sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi);
         -sin(theta), sin(phi) * cos(theta), cos(phi) * cos(theta)];
    elseif Flag == 'R'
        R = pose;
    end
end
    
        
 
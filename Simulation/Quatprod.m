function Qo = Quatprod(Q, Qr, flag)

% flag=0 表示单位四元数右乘向量或者四元数，要根据Qr的长度来确定计算方式
% 否则表示用四元数对向量做旋转，利用罗德里格斯公式可计算
    
    if flag ~= 0
        Qo = (eye(3) + 2 * R3so3(Q(2:4))^2 + 2 * Q(1) * R3so3(Q(2:4))) * Qr;
    else
        if size(Qr, 2) == 3
            Qo = [- Q(2:4)' * Qr;
                  cross(Q(2:4), Qr) + Q(1) * Qr];
        else
            Qo = [Q(1) * Qr(1) - Q(2:4)' * Qr(2:4);
                  cross(Q(2:4), Qr(2:4)) + Q(1) * Qr(2:4) + Qr(1) * Q(2:4)];
            Qo = Qo/norm(Qo);
        end
    end
    
end
    
    
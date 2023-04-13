function xhatNext = Smoother(u, xhat, t, para)
    global dt;
    if para(3) <= 0
        omegau = para(1);
    else
        omegau = TimeVaryingBuffer(t, para);
    end
    n = size(xhat, 1);
    A = zeros(n, n);
    A(1:n-1, 2:n) = eye(n-1);
%     A(2:n, 1:n-1) = omegau1d * diag(n-1:-1:1);
    B = zeros(n, 1);
    B(n) = omegau^n;
    A(n, n) = - n * omegau;
    for i = 2 : n
        A(n,n - i + 1) = A(n,n - i + 2) * (n - i + 1) / i * omegau;
    end
    
    xhatNext = xhat + dt * (A * xhat + B *(u));
end


% function xhatNext = Smoother(u, xhat, t, para)
%     global dt;
%     if para(3) <= 0
%         omegau = para(1);
%     else
%         omegau = TimeVaryingBuffer(t, para);
%     end
%     n = size(xhat, 1);
%     A = zeros(n, n);
%     A(1:n-1, 2:n) = eye(n-1);
% %     A(2:n, 1:n-1) = omegau1d * diag(n-1:-1:1);
%     L = zeros(n, 1);
%     L(1) = - n * omegau;
%     for i = 2 : n
%         L(i) = L(i-1) * (n - i + 1) / i * omegau;
%     end
%     
%     xhatNext = xhat + dt * (A * xhat + L *(xhat(1) - u));
% end

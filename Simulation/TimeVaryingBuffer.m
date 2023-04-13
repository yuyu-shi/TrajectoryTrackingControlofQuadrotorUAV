function omegau = TimeVaryingBuffer(t, para)
    omegauinfinity = para(1);
    omegau0 = para(2);
    tauk = para(3);
%     persistent omegauinfinity tauk omegau0;
%     if isempty(omegauinfinity)
%         omegauinfinity = 300;
%         tauk = 0.7;
%         if e == 0
%             omegau0 = 0;
%         else
%             omegau0 = 0.03;
%         end
%     end
    omegau = omegauinfinity / (1 + exp(-t/tauk) * (omegauinfinity/omegau0 - 1)); 

%     if t < tauk
%         omegau = omegau0 + (omegauinfinity - omegau0) * ((sin(pi * (t / tauk - 0.5)) + 1) / 2);
% %         omegau1d = (omegauinfinity - omegau0) / 2 * pi / tauk * cos(pi * (t / tauk - 0.5));
% %         omegau = omega0 / (1 + (omega0/omegauinfinity - 1)/2*(sin(pi * (t/tauk - 0.5)) + 1));
% %         omegau1d = 0;
% %         omegau = omegauinfinity / (1 + (omegauinfinity/omegau0 - 1) * (1 - t / tauk)^2);
% %         omegau1d = 0;
%     else 
%         omegau = omegauinfinity; 
% %         omegau1d = 0;
%     end
end

% function [pdhat, psidhat] = Smoother(u, tau)
%     global dt;
%     persistent uhat;
%     if isempty(uhat)
%         uhat = [0 0 0 0 0;
%                 0 0 0 0 0;
%                 0 0 0 0 0;
%                 0 0 0 0 0];
%     end
%     Au = [-1./tau.^5, -5./tau.^4, -10./tau.^3, -10./tau.^2, -5./tau];
%     Bu = 1./tau.^5;
%     uhat(:, 1:4) = uhat(:, 2:5) * dt + uhat(:, 1:4);
%     for i = 1 : 4
%         uhat(i, 5) = (Au(i, :) * uhat(i, :)' + Bu(i) * u(i)) * dt + uhat(i, 5);
%     end
%     pdhat = uhat(1:3, :);
%     psidhat = uhat(4, :);
% end
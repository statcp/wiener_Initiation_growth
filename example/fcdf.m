function FT = fcdf(t,w,mu0,sig02,a,b,sigb2)
%
% phi = mu0/sig0^2-2*(a-b*w)/sigb^2;
% psi = 1/sig0^2+4*b/sigb^2;
% 
% u0 = 2*a*w/sigb^2-mu0^2/(2*sig0^2)+phi^2/psi;
% u1 = ((1+b*t)*mu0-w+a*t)/(sigb^2*t+(1+b*t)^2*sig0^2)^.5;
% u2 = ((1-b*t)*phi-(w+a*t)*psi)/(sigb^2*t*psi^2+(1-b*t)^2*psi)^.5;
% 
% FT = normcdf(u1)+exp(u0)*normcdf(u2)/(sig0*psi^.5);



% mu0 = mu0G(j);sig0 = sig0G(j);a = aG(j);b = bG(j); sigb = sigbG(j);t=t0

phi = 2.*(b.*w-a)./sigb2+mu0./sig02;
psi = 4.*b./sigb2+1./sig02;

u1 = (a.*t+(1+b.*t).*mu0-w)./(sigb2.*t+(1+b.*t).^2.*sig02).^.5;
u2 = 2.*a.*w./sigb2-mu0.^2./(2.*sig02)+phi.^2./(2.*psi);
u3 = ((1-b.*t).*phi./psi-(w+a.*t))./(sigb2.*t+(1-b.*t).^2./psi).^.5;
FT = normcdf(u1)+exp(u2 + log(normcdf(u3)))./((sig02.*psi).^.5);


% phi = 2*(b*w-a)/sigb2+mu0/sig02;
% psi = 4*b/sigb2+1/sig02;
% 
% u1 = (a*t+(1+b*t)*mu0-w)/(sigb2*t+(1+b*t)^2*sig02)^.5;
% u2 = 2*a*w/sigb2-mu0^2/(2*sig02)+phi^2/(2*psi);
% u3 = ((1-b*t)*phi/psi-(w+a*t))/(sigb2*t+(1-b*t)^2/psi)^.5;
% FT = normcdf(u1)+exp(u2 + log(normcdf(u3)))/((sig02*psi)^.5);

% FT = normcdf(u1)+exp(u2)*normcdf(u3)/(sig0*psi^.5);


end
function rul = frul(t,w,mu0,sig02,a,b,sigb2)
% t: RUL, lk
% w: w-yk

eta = a+b*mu0;
sige = (sigb2*t+(b*t)^2*sig02)^.5;

u1 = (eta*t-w)/sige;
u2 = 2*w*eta/sigb2+2*(b*w)^2*sig02/sigb2^2;
u3 = -(2*b^2*sig02*w*t+sigb2*(eta*t+w))/(sigb2*sige);
rul = normcdf(u1)+exp(u2+log(normcdf(u3)));


%% t=Inf
% eta = a+b*mu0;
% sige = b*sig0;
% 
% u1 = eta/sige;
% u2 = 2*w*eta/sigb^2+2*(b*sig0*w)^2/sigb^4;
% u3 = -2*sige*w/sigb^2-eta/sige;
% rul = normcdf(u1)+exp(u2)*normcdf(u3);



end


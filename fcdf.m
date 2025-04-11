function FT = fcdf(t,w,mu0,sig02,a,b,sigb2)

phi = 2*(b*w-a)/sigb2+mu0/sig02;
psi = 4*b/sigb2+1/sig02;

u1 = (a*t+(1+b*t)*mu0-w)/(sigb2*t+(1+b*t)^2*sig02)^.5;
u2 = 2*a*w/sigb2-mu0^2/(2*sig02)+phi^2/(2*psi);
u3 = ((1-b*t)*phi/psi-(w+a*t))/(sigb2*t+(1-b*t)^2/psi)^.5;
FT = normcdf(u1)+exp(u2 + log(normcdf(u3)))/((sig02*psi)^.5);



end
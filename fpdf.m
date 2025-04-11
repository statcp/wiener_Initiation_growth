function ft = fpdf(t,w,mu0,sig02,a,b,sigb2)

psi0 = (1+b.*t).^2.*sig02+sigb2.*t;
phi0 = ((1+b.*t).*(w-a.*t).*sig02+mu0.*sigb2.*t)./psi0 ;
ft = (w-phi0)./(2.*pi.*t.^2.*psi0).^.5.*exp(-(w-a.*t-(1+b.*t).*mu0).^2./(2.*psi0));


end
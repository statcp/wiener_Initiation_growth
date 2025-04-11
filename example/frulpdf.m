function fr = frulpdf(l,wrul,mu0,sig02,a,b,sigb2)
% pdf of rul
% l:lk
%wrul:w-yk
vv = b.^2.*sig02.*l+sigb2;
fr = wrul.*exp(-(wrul-(a+b.*mu0).*l).^2./(2.*l.*vv))./(2.*pi.*l.^3.*vv).^.5;

end


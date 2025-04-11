function fr = frulpdf_m1(l,wrul,a,sigb)
% pdf of rul
% l:lk
%wrul:w-yk

fr = wrul.*exp(-(wrul-a.*l).^2./(2.*l.*sigb.^2))./(2.*pi.*l.^3.*sigb.^2).^.5;

end


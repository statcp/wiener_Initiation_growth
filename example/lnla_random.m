function lnl = lnla_random(x,tij,yij,n,mi)

% no random initial value with fixed drift a
% a follows N(mua,sa2)
% x(1) = sigma_a
% x(2) = sigma_B

% three parameters: siga2,sigB2,mua
% s2 = x(1).*tij+x(2);
% x(3) = sum(sum(yij./s2))/(n.*mi.*tij./s2);
% lnl = n.*mi.*log(s2)./2+sum(sum((yij-x(3).*tij).^2./(2.*tij.*s2)));

% two parameters: siga2,sigB2 
% s2 = x(1).*tij+x(2);
% mua = sum(sum(yij./s2))/(n.*mi.*tij./s2);
% lnl = n.*mi.*log(s2)./2+sum(sum((yij-mua.*tij).^2./(2.*tij.*s2)));

% two parameters: siga,sigB
s2 = x(1)^2.*tij^2+x(2)^2.*tij;
mua = sum(sum(yij./s2))/(n.*mi.*tij./s2);
lnl = n.*mi.*log(2*pi)+n.*mi.*log(s2)./2+sum(sum((yij-mua.*tij).^2./(2.*s2)));

end
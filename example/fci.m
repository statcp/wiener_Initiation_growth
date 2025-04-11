function [GG,Bayes] = fci(dat,B2,tij)

y0 = dat(1,:);
yij = diff(dat);

n = size(y0,2); % number of unit
mi = size(yij,1); % number of measurement
% tij = 500/24; %hours are transformed to days
n1n = (n-1)/n;
M = n*mi;
n2 = n^.5;
n12 = n/2;

% w = 5; % threshold

% MLEs, CIs of our method
mu0 = mean(y0);
sig02 = var(y0)*n1n; % sigma0^2
A = n*mi*tij;
B = sum(mi.*tij.*y0);
C = sum(mi.*tij.*y0.^2);
D = sum(sum(yij));
E = sum(sum(y0.*yij));
den = A*C-B^2;
a = (C*D-B*E)/den;
b = (A*E-B*D)/den;
sigb2 = sum(sum((yij-a.*tij-b.*y0.*tij).^2./tij))/M;
muiG = (a+b.*y0).*tij;
pam0 = [mu0, sig02, a, b, sigb2];


   % Bayes initial value
    mu0B(1) = mu0;
    sig02B(1) = sig02;
    aB(1) = a;
    bB(1) = b;
    sigb2B(1) = sigb2;

    % GPQ(1)    
    U1 = chi2rnd(n-1,1);
    U2 = normrnd(0,1);
    sig02G(1) = n*sig02/U1;
    mu0G(1) =  mu0-U2*sig02G(1)^.5/n2;
    U3 = chi2rnd(M-2);
    sigb2G(1) = M*sigb2/U3;
    U4 = normrnd(0,1);
    U5 = normrnd(0,1);    
    bG(1) = b-U5*(sigb2G(1)*A/den)^.5;
    aG(1) = a+(B/A)*(b-bG(1))-U4*(sigb2G(1)/A)^.5;    
    
    
for j = 2:B2
%     j
     %% objective Bayesian
        mu0B(j) = normrnd(mu0,sig02B(j-1)^.5/n2);
        bays_ss = sum((y0-mu0B(j)).^2)/2;
        sig02B(j) = 1/gamrnd(n12,1/bays_ss);
        aB(j) = normrnd((D-bB(j-1)*B)/A, (sigb2B(j-1)/A)^.5);
        bB(j) = normrnd((E-aB(j)*B)/C, (sigb2B(j-1)/C)^.5);
        v3 = sum(sum((yij - aB(j)*tij-bB(j)*tij*y0).^2./tij))/2;
        sigb2B(j) = 1/gamrnd(M/2+1,1/v3);     
%         reB(j) = 1-fcdf(t0,w,mu0B(j),sig02B(j),aB(j),bB(j),sigb2B(j));
%         if 4*bB(j)/sigb2B(j)+1/sig02B(j) < 0
%             reB(j) = 1-integral(@(x) fpdf(x,w,mu0B(j),sig02B(j),aB(j),bB(j),sigb2B(j)),0,t0);
%         else reB(j) = 1-fcdf(t0,w,mu0B(j),sig02B(j),aB(j),bB(j),sigb2B(j));
%         end

        %% GPQ
        U1 = chi2rnd(n-1,1);
        U2 = normrnd(0,1);
        sig02G(j) = n*sig02/U1;
        mu0G(j) =  mu0-U2*sig02G(j)^.5/n2;
        U3 = chi2rnd(M-2);
        sigb2G(j) = M*sigb2/U3;
        U4 = normrnd(0,1);
        U5 = normrnd(0,1);
        bG(j) = b-U5*(sigb2G(j)*A/den)^.5;
        aG(j) = a+(B/A)*(b-bG(j))-U4*(sigb2G(j)/A)^.5; 
    
end
mu0B2 = mu0B((0.5*B2+1):end);
sig0B2 = sig02B((0.5*B2+1):end);
aB2 = aB((0.5*B2+1):end);
bB2 = bB((0.5*B2+1):end);
sigbB2 = sigb2B((0.5*B2+1):end);
    
gpq = [mu0G;sig02G;aG;bG;sigb2G];
bayes = [mu0B2;sig0B2;aB2;bB2;sigbB2];
 
GG = mean(gpq,2);
Bayes = mean(bayes,2);


end


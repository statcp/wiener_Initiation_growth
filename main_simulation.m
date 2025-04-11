% Wiener ramdon initial
clc;
clear;

%% Parameter setting
n = 10; % number of unit
mi = 10; % number of measurement
tij = 1;
mu0 = 3; sig02 = 0.25; a = 3; b = 1; sigb2 = .25; w = 45; % threshold

%% Useful variables
A = n*mi*tij;
M = n*mi;
n2 = n^.5;
n12 = n/2;
n1n = (n-1)/n;
options=optimset('Display','off');
re = 0.8;
t0 = fsolve(@(x) fcdf(x,w,mu0,sig02,a,b,sigb2)-1+re,6,options)

%% RUL 
pt = [1/3 1/2 2/3];
yk = pt.*w;
wrul = w-yk;
prul = size(pt,2);

%% Output setting
B1 = 10000;
B2 = 10000;
quan = [0.025,0.975].*100;
true0 = [mu0,sig02,a,b,sigb2,re];
true02 = repmat(true0,1,2);
CP = zeros(size(true02,2)+prul*2,1);Length=CP;


[mu0,sig02,a,b,sigb2,w] 
[n,mi]

for i = 1:B1
    i
    tic;

    %% generate samples
    y0 = normrnd(mu0,sig02^.5,1,n);
    mut = (a+b.*y0).*tij;
    sigt = (sigb2.*tij)^.5;
    for ii = 1:n
        yij(1:mi,ii) = normrnd(mut(ii),sigt,mi,1);
    end
    
    %% MLE
    m0mle = mean(y0);
    s02mle = var(y0)*n1n; % sigma0^2
    D = sum(sum(yij));
    B = sum(mi.*tij.*y0);
    C = sum(mi.*tij.*y0.^2);
    E = sum(sum(y0.*yij));
    den = A*C-B^2;
    amle = (C*D-B*E)/den;
    bmle = (A*E-B*D)/den;
    sb2mle = sum(sum((yij-amle.*tij-bmle.*y0.*tij).^2./tij))/M;
    
     %% Bayes initial value
     mu0B(1) = m0mle;
     sig02B(1) = s02mle;
     aB(1) = amle;
     bB(1) = bmle;
     sigb2B(1) = sb2mle;
     if 4*bB(1)/sigb2B(1)+1/sig02B(1) < 0
         reB(1) = 1-integral(@(x) fpdf(x,w,mu0B(1),sig02B(1),aB(1),bB(1),sigb2B(1)),0,t0);
     else reB(1) = 1-fcdf(t0,w,mu0B(1),sig02B(1),aB(1),bB(1),sigb2B(1));
     end
     % GPQ(1)
     U1 = chi2rnd(n-1,1);
     U2 = normrnd(0,1);
     sig02G(1) = n*s02mle/U1;
     mu0G(1) =  m0mle-U2*sig02G(1)^.5/n2;
     U3 = chi2rnd(M-2);
     sigb2G(1) = M*sb2mle/U3;
     U4 = normrnd(0,1);
     U5 = normrnd(0,1);
     bG(1) = bmle-U5*(sigb2G(1)*A/den)^.5;
     aG(1) = amle+(B/A)*(bmle-bG(1))-U4*(sigb2G(1)/A)^.5;
     if 4*bG(1)/sigb2G(1)+1/sig02G(1) < 0
        reG(1) = 1-integral(@(x) fpdf(x,w,mu0G(1),sig02G(1),aG(1),bG(1),sigb2G(1)),0,t0);
    else reG(1) = 1-fcdf(t0,w,mu0G(1),sig02G(1),aG(1),bG(1),sigb2G(1));
    end
     %% RUL    
     R0 = unifrnd(0,1,1);
     R1 = unifrnd(0,1,1);
     for ii = 1:prul
         rul(ii) = fbrul(R0,wrul(ii),mu0,sig02,a,b,sigb2); % true values of RULs
         rulB(1,ii) = fbrul(R1,wrul(ii),mu0B(1),sig02B(1),aB(1),bB(1),sigb2B(1));
         rulG(1,ii) = fbrul(R1,wrul(ii),mu0G(1),sig02G(1),aG(1),bG(1),sigb2G(1));
     end
     
    %%
    for j = 2:B2
       %% objective Bayesian
        mu0B(j) = normrnd(m0mle,sig02B(j-1)^.5/n2);
        bays_ss = sum((y0-mu0B(j)).^2)/2;
        sig02B(j) = 1/gamrnd(n12,1/bays_ss);
        aB(j) = normrnd((D-bB(j-1)*B)/A, (sigb2B(j-1)/A)^.5);
        bB(j) = normrnd((E-aB(j)*B)/C, (sigb2B(j-1)/C)^.5);
        v3 = sum(sum((yij - aB(j)*tij-bB(j)*tij*y0).^2./tij))/2;
        sigb2B(j) = 1/gamrnd(M/2+1,1/v3);
        if 4*bB(j)/sigb2B(j)+1/sig02B(j) < 0
            reB(j) = 1-integral(@(x) fpdf(x,w,mu0B(j),sig02B(j),aB(j),bB(j),sigb2B(j)),0,t0);
        else reB(j) = 1-fcdf(t0,w,mu0B(j),sig02B(j),aB(j),bB(j),sigb2B(j));
        end
        %% GPQ
        U1 = chi2rnd(n-1,1);
        U2 = normrnd(0,1);
        sig02G(j) = n*s02mle/U1;
        mu0G(j) =  m0mle-U2*sig02G(j)^.5/n2;
        U3 = chi2rnd(M-2);
        sigb2G(j) = M*sb2mle/U3;
        U4 = normrnd(0,1);
        U5 = normrnd(0,1);
        bG(j) = bmle-U5*(sigb2G(j)*A/den)^.5;
        aG(j) = amle+(B/A)*(bmle-bG(j))-U4*(sigb2G(j)/A)^.5;
        if 4*bG(j)/sigb2G(j)+1/sig02G(j) < 0
            reG(j) = 1-integral(@(x) fpdf(x,w,mu0G(j),sig02G(j),aG(j),bG(j),sigb2G(j)),0,t0);
        else reG(j) = 1-fcdf(t0,w,mu0G(j),sig02G(j),aG(j),bG(j),sigb2G(j));
        end
        %% RUL
        R2 = unifrnd(0,1,1);
        for ii = 1:prul
            rulG(j,ii) = fbrul(R2,wrul(ii),mu0G(j),sig02G(j),aG(j),bG(j),sigb2G(j));
            rulB(j,ii) = fbrul(R2,wrul(ii),mu0B(j),sig02B(j),aB(j),bB(j),sigb2B(j));
        end
    end

      
    % quantile of objective Bayesian
    mu0B2 = mu0B((0.5*B2+1):end);
    sig0B2 = sig02B((0.5*B2+1):end);
    aB2 = aB((0.5*B2+1):end);
    bB2 = bB((0.5*B2+1):end);
    sigbB2 = sigb2B((0.5*B2+1):end);    
    reB2 = reB((0.5*B2+1):end);
    qmu0B = prctile(mu0B2,quan); 
    qsig02B = prctile(sig0B2,quan);
    qaB = prctile(aB2,quan);
    qbB = prctile(bB2,quan);
    qsigb2B = prctile(sigbB2,quan);      
    qreB = prctile(reB2,quan);
    qbayes = [qmu0B;qsig02B;qaB;qbB; qsigb2B;qreB];
    rulB2 = rulB((0.5*B2+1):end,:);
    qrulB = prctile(rulB2,quan)'; 
    
    % quantile of GPQ
    qmu0G = prctile(mu0G,quan); 
    qsig02G = prctile(sig02G,quan);
    qaG = prctile(aG,quan);
    qbG = prctile(bG,quan);
    qsigb2G = prctile(sigb2G,quan);
    qreG = prctile(reG,quan);
    qgpq = [qmu0G; qsig02G; qaG;qbG; qsigb2G;qreG];
    qrulG = prctile(rulG,quan)';  
    
    sort_all = [qbayes;qgpq;qrulB;qrulG];
    true = [true02,rul,rul];


    %% Coverage probabilities
    for row = 1:size(sort_all,1)
        if (sort_all(row,1)<= true(row))  &&  (true(row)<=sort_all(row,2))
            CP(row,1)=CP(row,1)+1;
        end
        Length(row,1)=Length(row,1)+sort_all(row,2)-sort_all(row,1);  
    end 


toc;

end


CP2 = CP/B1;
Length2 = Length/B1;

sn = size(true0,2);
bayescp = CP2(1:sn,:);
gpqcp = CP2(1+sn:2*sn,:);
bayesrulcp = CP2(1+2*sn:3+2*sn,:);
gpqrulcp = CP2(4+2*sn:end,:);

bayeslen = Length2(1:sn,:);
gpqlen = Length2(1+sn:2*sn,:);
bayesrullen = Length2(1+2*sn:3+2*sn,:);
gpqrullen = Length2(4+2*sn:end,:);


[mu0,sig02,a,b,sigb2,re] 
[n,mi]
bayescp,bayeslen
gpqcp,gpqlen
bayesrulcp,bayesrullen
gpqrulcp,gpqrullen



% save("a4-10-20")
% t2 = datestr(clock)
% clc;clear





















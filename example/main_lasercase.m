% Wiener ramdon initial
clc;
clear;
% M0: our model
% M1: fixed effect Wiener
% M2: independent
% M3: random effect
  tic;
%% Laser Data
DD = [
0.4741	0.7078	0.7074	0.3634	0.2707	0.3584	0.3637	0.458	0.5098	0.4136	0.4435	0.3864	0.3034	0.438	0.5056
0.9255	1.2175	1.1651	0.6161	0.6144	1.3922	0.9237	1.0678	0.9308	1.488	0.9978	0.7988	0.7448	0.6979	0.8349
2.1147	1.8955	1.7253	1.3597	1.1093	1.9541	1.2098	1.4163	1.5713	2.381	1.5738	1.3521	1.5225	1.0515	1.2864
2.7168	2.2954	1.9888	1.9525	1.773	2.8617	1.4591	1.7659	1.9576	2.995	1.9584	1.741	1.8548	1.3526	1.5191
3.511	2.865	2.5325	2.2995	2.0625	3.4552	1.9323	2.1061	2.5871	3.835	2.5093	2.979	2.3897	1.7998	1.9148
4.3415	3.7454	2.9695	2.9468	2.5835	3.8117	2.3862	2.3969	3.2892	4.501	2.8436	3.5871	2.9473	2.5514	2.2697
4.9076	4.4192	3.2977	3.3871	2.9863	4.533	2.6799	2.78	3.6075	5.251	3.4682	4.0315	3.5104	2.8274	2.7763
5.4782	4.9894	3.9354	3.7863	3.3765	5.3541	2.9397	3.0187	4.1148	6.256	4.007	4.4393	3.9207	3.3908	3.4153
5.9925	5.5148	4.1613	4.1116	4.0491	5.9217	3.4183	3.289	4.5964	7.051	4.5124	4.7906	5.0316	3.7159	3.7814
6.7224	6.0668	4.4459	4.5013	4.626	6.7078	4.0857	3.7542	4.913	7.803	4.8018	5.2205	5.4708	4.0906	4.1053
7.1303	6.6385	4.889	4.7209	5.237	7.7027	4.5843	4.1595	5.3408	8.321	5.2042	5.4848	5.8362	4.8296	4.3754
8.0006	7.1615	5.2696	4.9766	5.6191	8.6053	4.8402	4.7575	5.8418	8.93	5.6556	5.9639	6.4962	5.4145	4.6295
8.9193	7.7778	5.6913	5.2812	6.0449	9.1477	5.1145	5.1597	6.3961	9.554	6.1951	6.2276	6.9372	5.759	5.3831
9.494	8.4242	6.0216	5.615	6.3167	9.951	5.5691	5.4591	6.8441	10.45	6.5405	6.9877	7.3883	6.1419	5.8391
9.8675	8.9092	6.4485	5.9503	7.0952	10.4857	6.114	5.8085	7.2037	11.28	6.9626	7.3671	7.8472	6.5086	6.1618
10.9446	9.2834	6.8849	6.1438	7.5941	11.0096	7.1723	6.245	7.8848	12.21	7.4238	7.8841	8.087	6.8823	6.6247
 ];
% 15 unit are tested at 15 time points, row is units, col is time
w = 10; % threshold
T1 = [250 500 750	1000 1250 1500 1750	2000 2250 2500 2750	3000 3250 3500 3750	4000];
burn_in = 500;
location_burnin = find(T1 == burn_in );
dat0 = DD(1:location_burnin-1,:);
T0 = T1(1:location_burnin-1);
dat = DD(location_burnin:end,:); %mi*n
TT = T1(location_burnin:end)
tij = 250;

%% Figure 10 in the case study
% Degradation paths for the GaAs Laser data with 500 hours burn-in period.
figure;
hold on;
legendHandles = []; 
legendInfo = {};  
numColors = size(DD, 2); 
hues = linspace(0, 1, numColors+1);
hues(end) = [];  
colors = hsv2rgb([hues' ones(numColors, 1) ones(numColors, 1)]); 
lightColors = colors + 0.7 * (1 - colors); 
darkColors = 1 * colors; 

markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', '*', '+', 'x', 'h', '.', '|', '_', 'o'};
% T1 = T1./24;
for i = 1:size(DD, 2)
    markerIndex = mod(i-1, length(markers)) + 1;
    marker = markers{markerIndex};
    plot(T1, DD(:, i), [marker '-'], 'Color', lightColors(i, :), 'LineWidth', 1, 'MarkerSize', 4);
    p2 = plot(T1(2:end), DD(2:end, i), [marker '-'], 'Color', darkColors(i, :), 'LineWidth', 1, 'MarkerSize', 4);
    if isvalid(p2)
        legendHandles = [legendHandles, p2];
        legendInfo{end+1} = ['# ' num2str(i)]; 
    end
end
hold on;
box on;
xline(500,'r--','linewidth',2)
text(600,4.5,'End of burn in: 500 hr','Rotation',90,'Color','red')
xlabel('Hours')
yline(10,'linewidth',1)
text(600,11,'Threshold')
ylabel('Percent Increase in Operating Current')



%% Table 1 and Table 2 in the case study
% Intervals of M0
y0 = dat(1,:);
yij = diff(dat);
%  Useful variables
n = size(y0,2); % number of unit
mi = size(yij,1); % number of measurement
n1n = (n-1)/n;
M = n*mi;
n2 = n^.5;
n12 = n/2;
tij = tij./24;
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
options=optimset('Display','off');
% RUL
pt = [1/3 1/2 2/3];
yk = pt.*w;
wrul = w-yk;
prul = size(pt,2);

% Output setting
B2 = 10000;
quan = [0.025,0.975].*100;

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

% RUL
R1 = unifrnd(0,1,1);
for ii = 1:prul
    rulB(1,ii) = fbrul(R1,wrul(ii),mu0B(1),sig02B(1),aB(1),bB(1),sigb2B(1));
    rulG(1,ii) = fbrul(R1,wrul(ii),mu0G(1),sig02G(1),aG(1),bG(1),sigb2G(1));
end
    
for j = 2:B2
 
        % objective Bayesian
        mu0B(j) = normrnd(mu0,sig02B(j-1)^.5/n2);
        bays_ss = sum((y0-mu0B(j)).^2)/2;
        sig02B(j) = 1/gamrnd(n12,1/bays_ss);
        aB(j) = normrnd((D-bB(j-1)*B)/A, (sigb2B(j-1)/A)^.5);
        bB(j) = normrnd((E-aB(j)*B)/C, (sigb2B(j-1)/C)^.5);
        v3 = sum(sum((yij - aB(j)*tij-bB(j)*tij*y0).^2./tij))/2;
        sigb2B(j) = 1/gamrnd(M/2+1,1/v3);     

        % GPQ
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

        % RUL
        R1 = unifrnd(0,1,1);
        for ii = 1:prul
            rulG(j,ii) = fbrul(R1,wrul(ii),mu0G(j),sig02G(j),aG(j),bG(j),sigb2G(j));
            rulB(j,ii) = fbrul(R1,wrul(ii),mu0B(j),sig02B(j),aB(j),bB(j),sigb2B(j));
        end
end

    % quantile of objective Bayes
    mu0B2 = mu0B((0.5*B2+1):end);
    sig0B2 = sig02B((0.5*B2+1):end);
    aB2 = aB((0.5*B2+1):end);
    bB2 = bB((0.5*B2+1):end);
    sigbB2 = sigb2B((0.5*B2+1):end);    
    qmu0B = prctile(mu0B2,quan); 
    qsig02B = prctile(sig0B2,quan);
    qaB = prctile(aB2,quan);
    qbB = prctile(bB2,quan);
    qsigb2B = prctile(sigbB2,quan);      
    
     % quantile of GPQ
    qmu0G = prctile(mu0G,quan); 
    qsig02G = prctile(sig02G,quan);
    qaG = prctile(aG,quan);
    qbG = prctile(bG,quan);
    qsigb2G = prctile(sigb2G,quan);

    rulB2 = rulB((0.5*B2+1):end,:);
    qrulB = prctile(rulB2,quan)';  
    qrulG = prctile(rulG,quan)';  

% Table 1 in the case study
mu0
[qmu0B; qmu0G]
sig02.*100
[qsig02B;qsig02G].*100
a.*100
[qaB;qaG].*100
b.*100
[qbB;qbG].*100
sigb2.*1000
[qsigb2B;qsigb2G].*1000

% Table 2 in the case study
qrulB, qrulG


%% Figure 11 in the case study
% RUL of GPQ and Bayes
TT2 = TT./24;
Trul = (2500/24):tij:(4000/24);
low_t1 = find(TT2 == Trul(1));
up_t1 = find(TT2 == Trul(end));
ran_t1 = round(max(diff(Trul)/tij));
unit10 = dat(:,1);
rul_true = 4000/24-Trul %unit10
unit = unit10;
yk = unit((low_t1-1):ran_t1:(up_t1-1)); %current at 750,1250,1750,2250,2750
wrul = w-yk;
B3 = 5000;
clear GG Bayes
for kk = 1:size(wrul,1)
[GG(:,kk),Bayes(:,kk)] = fci(dat(1:(low_t1-1+(kk-1)*ran_t1),:),B3,tij);
end

% clear R1 bayesR1 tt
tt = 0.00001:1:(4000/24);
% GPQ and BP based on PDF of RUL
for pp = 1:size(wrul,1)
    for i = 1:size(tt,2)
        t = tt(i);        
        R1(i,pp) = frulpdf(t,wrul(pp),GG(1,pp),GG(2,pp),GG(3,pp),GG(4,pp),GG(5,pp));
        bayesR1(i,pp) = frulpdf(t,wrul(pp),Bayes(1,pp),Bayes(2,pp),Bayes(3,pp),Bayes(4,pp),Bayes(5,pp));
    end
end
% plot(R1)
[max_R1,meangpq]=max(R1,[],1)
[max_bayesR1,meanbayes]=max(bayesR1,[],1)

si1 = size(tt,1);
si2 = size(tt,2);
onett = ones(si2,si1);
TT3 = Trul.*onett;

p1 = plot3(tt,TT3(:,1),R1(:,1),'LineWidth',1,'Color','r')
hold on
p2 = plot3(tt,TT3,R1,'LineWidth',1,'Color','r')
p3 = plot3(tt,TT3(:,1),bayesR1(:,1),'LineWidth',1,'Color','b')
p4 = plot3(tt,TT3,bayesR1,'LineWidth',1,'Color','b')
p5 = plot(meangpq,Trul, '*r','LineWidth',1) 
p6 = plot(meanbayes,Trul, '^b','LineWidth',1) 
p7 = plot(rul_true,Trul, '--k','LineWidth',1) % true RUL
hold off
legend([p1 p3 p7],{'Generalized method','Bayesian method','True RUL'},'Location','north' )
grid on
xlabel('RUL')
ylabel('Days')
zlabel('PDF of RUL')
xlim([0,round(4000/24)])
ylim([2500/24,round(4000/24)])



%% Table 3 in the case study
% Log-likelihood and AIC values of M0,M1,M2,M3

% M0
pam0 = [mu0, sig02, a, b, sigb2];
muiG = (a+b.*y0).*tij;
lnl_m0 = -(n/2)*log(2*pi*sig02)-sum((y0-mu0).^2./(2*sig02))-.5*M*log(2*pi*tij*sigb2)-sum(sum((yij-muiG).^2./(2*sigb2*tij)));
AIC_m0 = 2*5-2*lnl_m0;

% M1：y0 = 0, Y(t) = at+sigb B(t)
yij_m1 = yij;
a_m1 = sum(sum(yij_m1))/(M*tij);
sigb2_m1 = sum(sum((yij_m1-a_m1*tij).^2./tij))/M;
pam1 = [a_m1,sigb2_m1];
% a = sum(range(dat))/(n*(4000-250))
lnl_m1 = -(M/2)*log(2*pi*tij*sigb2_m1)-sum(sum((yij_m1-a_m1*tij).^2./(2*sigb2_m1*tij)));
AIC_m1 = 2*2-2*lnl_m1;

% M2: b = 0,y0~N(mu0,s02), Y(t) = y0+at+sigb B(t)
mu0_m2 = mu0;
sig02_m2 = sig02;
a_m2 = sum(sum(yij))/(M*tij);
sigb2_m2 = sum(sum((yij-a_m2.*tij).^2./tij))/M;
pam2 = [mu0_m2,sig02_m2,a_m2,sigb2_m2];
lnl_m2 = -(n/2)*log(2*pi*sig02_m2)-sum((y0-mu0_m2).^2./(2*sig02_m2))-.5*M*log(2*pi*tij*sigb2_m2)-sum(sum((yij-a_m2*tij).^2./(2*sigb2_m2*tij)));
AIC_m2 = 2*4-2*lnl_m2;

% M3:Wiener with ramdom effect, Y(t) = at+sigb B(t), a~N(mua,siga2)
x0 = [sig02,sigb2];
sig_m3 = fminsearch(@(x) lnla_random(x,tij,yij,n,mi),x0,options);
sa2 = sig_m3(1).^2;
sigb2_m3 = sig_m3(2).^2;
s2_m3 = sa2.*tij+sigb2_m3;
mua = sum(sum(yij./s2_m3))/(n.*mi.*tij./s2_m3);
pam3 = [mua,sig_m3(1),sig_m3(2)];
lnl_m3 = -(M/2)*log(2*pi*tij*s2_m3)-sum(sum((yij-mua.*tij).^2./(2.*tij.*s2_m3)));
AIC_m3 = 2*3-2*lnl_m3;


% Table 3
lnl = [lnl_m0, lnl_m1, lnl_m2, lnl_m3]
AIC = [AIC_m0, AIC_m1, AIC_m2, AIC_m3]


%% Figure 12 in the case study
% RUL of M0, M1
Trul = 1000:500:4000;
Trul = Trul./24;
low_t1 = find(TT2== Trul(1));
up_t1 = find(TT2 == Trul(end));
rul_true = 4000/24-Trul; %unit1
unit = dat(:,1);
yk = unit((low_t1-1):2:(up_t1-1)); %current at 750,1250,1750,2250,2750
wrul = w-yk;

clear tt rulm0 rulm1
tt = 0:1:300;
% GPQ and BP based on all data
for pp = 1:size(wrul,1)
    for i = 1:size(tt,2)
        t=tt(i);
        rulm0(i,pp) = frulpdf(t,wrul(pp),mu0,sig02,a,b,sigb2);
        rulm1(i,pp) = frulpdf_m1(t,wrul(pp),a_m1,sigb2_m1^.5);
    end
end
[max_rulm0,meanm0]=max(rulm0,[],1)
[max_rulm1,meanm1]=max(rulm1,[],1)


si1 = size(tt,1);
si2 = size(tt,2);
onett = ones(si2,si1);
TT4 = Trul.*onett;

p01 = plot3(tt,TT4(:,1),rulm0(:,1),'LineWidth',1,'Color','r');
hold on
p02 = plot3(tt,TT4,rulm0,'LineWidth',1,'Color','r');
p03 = plot3(tt,TT4(:,1),rulm1(:,1),'LineWidth',1,'Color','b');
p04 = plot3(tt,TT4,rulm1,'LineWidth',1,'Color','b');
p05 = plot(meanm0,Trul, '*r','LineWidth',1) ;
p06 = plot(meanm1,Trul, '^b','LineWidth',1) ;
% legend([p01 p03],{'M0','M1'},'Location','north' )
p07 = plot(rul_true,Trul, '--k','LineWidth',1) % true RUL
hold off
legend([p01 p03 p07],{'M0','M1','True RUL'},'Location','north' )
grid on
xlabel('RUL')
ylabel('Days')
zlabel('PDF of RUL')





toc;



clear all
hold off
timesteps = 36;
lambda = 1;
X = linspace(0.0,10,500);
b = exp(-X)./(2*pi*besseli(0,X));
chifactor = pi;
fidi = fopen('bprw_full_CI_CR_julien.txt');
t = textscan(fidi, '%f%f%f%f%f%f%f%f', 'Delimiter',',');
kappabmax = 10;
kappa1max = 0.2;
kappamax = 10;
alphamin = -20; 
alphamax = 20;
fclose(fidi);
alpha=t{1};
kappa=t{2};
kappab=t{3};
kappa1=t{4};
CR=t{5};
CI=t{6};
CR_analytic=t{7};
CI_analytic=t{8};

%x_analytic,rsqr_analytic,distavg,distsqravg,xavg,xsqavg,xrratio,yavg
fidi = fopen('bprw_CI_CR_julien.txt');
y = textscan(fidi, '%f%f%f%f%f%f%f%f', 'Delimiter',',');
fclose(fidi);
x_analytic=y{1};
rsqr_analytic=y{2};
distavg=y{3};
distsqravg=y{4};
xavg=y{5};
xsqavg=y{6};
xrratio=y{7};
yavg=y{8};

CR_exp1 = 0.546347;
CR_exp2 = 0.523874;
CR_exp3 = 0.507164;
CR_exp4 = 0.472901;

dCR_exp1 = 0.043344;
dCR_exp2 = 0.047813;
dCR_exp3 = 0.039485;
dCR_exp4 = 0.054968;

dCI_exp1 = 0.107776;
dCI_exp2 = 0.035738;
dCI_exp3 = 0.095798;
dCI_exp4 = 0.108965;

CI_exp1 = -0.09441;
CI_exp2 = -0.01209;
CI_exp3 = 0.309924;
CI_exp4 = 0.342912;


dd2 = '2018_04_02_tau_100_sec_alpha_2_lambda_10_r_1_by_tau_100pts_N_10000_z_100_c_codes_cpm6s/dat/';

%param = load([dd2 'cpm5g.param.dat']);
g21(1)=0.01;
g21(2)=1;
g21(3)=5;
g21(4)=50;


    for j = 1:4
        v10 = load([dd2 'cpm5g_g' num2str(j-1) '.v.dat']);
        CI10 = load([dd2 'cpm5g_g' num2str(j-1) '.CI.dat']);
        CR10 = load([dd2 'cpm5g_g' num2str(j-1) '.CR.dat']);

        v10 = v10*60*60; % um/h

        Z = length(v10);
        vbar2(j) = mean(v10);
        dv2(j) = std(v10)/sqrt(Z);
        verr2(j)=log10(vbar2(j)+dv2(j))-log10(vbar2(j));
        CIbar2(j) = mean(CI10);
        dCI2(j) = std(CI10)/sqrt(Z);
        CRbar2(j) = mean(CR10);
        dCR2(j) = std(CR10)/sqrt(Z);
    end

X1 = linspace(0.0,100,50000);
hold off;
b1 = exp(-X1)./(2*pi*besseli(0,X1));
%X = linspace(0.0,2,1000);
%psi = besseli(1,X)./besseli(0,X);
%b = 0.16;
psi = besseli(1,X1)./besseli(0,X1);
chi = chifactor*b1;
%t = 36;
%lambda = 1;
N = lambda*timesteps;
g1 = timesteps*(1-((2*chi.^2)./(1-psi).^2));
g2 = -timesteps*(chi.^2).*exp(-N*(1-psi))./(1-psi).^2;
g3 = (2*(chi.^2)./(1-psi).^2-1)./(lambda*(1-psi)).*(1-exp(-N*(1-psi)));
g4 = 0.5*(chi.^2)*N*timesteps./(1-psi);
g5 = ((chi.^2).*(1-exp(-N*(1-psi)).^2)./(lambda*(1-psi).^3));
g = g1+g2+g3+g4+g5;

CI1 = chi.*sqrt(0.5*lambda./(1-psi)).*(timesteps-(1-exp(-N*(1-psi)))/(lambda.*(1-psi)))./sqrt(g);
CR1 = (1/timesteps)*sqrt(2*g./(lambda*(1-psi)));


X2 = 0;
hold off;
b2 = linspace(0.0001,0.5/pi,100000);
%X = linspace(0.0,2,1000);
%psi = besseli(1,X)./besseli(0,X);
%b = 0.16;
psi = besseli(1,X2)/besseli(0,X2);
chi = chifactor*b2;
%t = 36;
%lambda = 1;
N = lambda*timesteps;
g1 = timesteps*(1-((2*chi.^2)/(1-psi).^2));
g2 = -timesteps*(chi.^2)*exp(-N*(1-psi))/(1-psi).^2;
g3 = (2*(chi.^2)/(1-psi)^2-1)/(lambda*(1-psi))*(1-exp(-N*(1-psi)));
g4 = 0.5*(chi.^2)*N*timesteps/(1-psi);
g5 = ((chi.^2)*(1-exp(-N*(1-psi))^2)/(lambda*(1-psi)^3));
g = g1+g2+g3+g4+g5;

CI2 = chi.*sqrt(0.5*lambda./(1-psi)).*(timesteps-(1-exp(-N*(1-psi)))/(lambda.*(1-psi)))./sqrt(g);
CR2 = (1/timesteps)*sqrt(2*g./(lambda*(1-psi)));

X3 = linspace(0.0,100,100000);
hold off;
b3 = 0;
%X = linspace(0.0,2,1000);
%psi = besseli(1,X)./besseli(0,X);
%b = 0.16;
psi = besseli(1,X3)./besseli(0,X3);
chi = chifactor*b3;
%t = 36;
%lambda = 1;
N = lambda*timesteps;
g1 = timesteps*(1-((2*chi^2)./(1-psi).^2));
g2 = -timesteps*(chi^2).*exp(-N*(1-psi))./(1-psi).^2;
g3 = (2*(chi^2)./(1-psi).^2-1)./(lambda*(1-psi)).*(1-exp(-N*(1-psi)));
g4 = 0.5*(chi^2)*N*timesteps./(1-psi);
g5 = ((chi^2)*(1-exp(-N*(1-psi)).^2)./(lambda*(1-psi).^3));
g = g1+g2+g3+g4+g5;

CI3 = chi.*sqrt(0.5*lambda./(1-psi)).*(timesteps-(1-exp(-N*(1-psi)))/(lambda.*(1-psi)))./sqrt(g);
CR3 = (1/timesteps)*sqrt(2*g./(lambda*(1-psi)));



hold off;


%t = 36;
%lambda = 1;
N = lambda*timesteps;
CR4 = linspace(sqrt(2/(N+2)),sqrt(0.25+1/N),1000);

%CI4 = sqrt(1-(2/N)*(1./CR4.^2));
CI4 = sqrt(1-(2/N)*((1-CR4.^2)./CR4.^2));

N = lambda*timesteps;
CR5 = linspace(0,1,1000);

CI5 = 2*N*(CR5.^2).*exp(-4*N*CR5.^2)*sqrt(pi*N/16);

N = lambda*timesteps;
CR6 = linspace(0.3,1,1000);
p1 = 9*CR6./(32*(1-CR6));
bracket1 = 0.5*N;
bracket2 = sqrt(0.5*pi*p1);
bracket3 = exp(-2*p1).*(1-13./(32*p1)-(9./(256*p1.^2)));
CI6 = bracket1.*bracket2.*bracket3;

N = lambda*timesteps;
p2 = linspace(1,100,1000);
bracket1 = 16*p2.^2;
bracket2 = 16*p2.^2+9*p2+(81/64);
CRtest = sqrt(bracket1./bracket2);

CR7 = linspace(0.3,1,1000);
alpha = 0.125*N^2+0.5*N-0.25;
beta = -(3/16)*N^2+0.5*N-13/8;
gamma = N^2/32+N/8;
pa = beta/2+gamma;
pb = alpha/2+beta;
pc = alpha-(N^2)*(CR7.^2)/2;
p7pos = (-pb+sqrt(pb^2-4*pa*pc))/(2*pa);
p7neg = (-pb-sqrt(pb^2-4*pa*pc))/(2*pa);
CI7pos = ((1-1/N)-(1-2/N)*p7pos/2+(p7pos.^2)/4)./(2*CR7);
CI7neg = ((1-1/N)-(1-2/N)*p7neg/2+(p7neg.^2)/4)./(2*CR7);



CR8 = linspace(0,1,1000);
pa = (1/6)*(N./(1-CR8.^2)-(27/16));
CI8 = (1./CR8).*sqrt(0.5*pi*pa).*exp(-2*pa).*(0.5*N-(N^2/12+1/8)./pa-(5*N^2/384+9*N/256)./(pa.^2));

N = lambda*timesteps;
CR9 = linspace(0.3,1,1000);
alpha = (N^2+4*N-2)/(4*N^2);
beta = 2*(-1/N^2+3/(4*N)-1/8);
gamma = (N^2+8*N-24)/(16*N^2);
p9 = (-beta+sqrt(beta^2+4*gamma*(alpha-CR9.^2)))/(2*gamma);
M9 = N*(1-0.5*p9);
z9 = 0.5-0.25*p9;
CI9 = (z9./CR9).*(M9-1)./M9;

N = lambda*timesteps;
CR10 = linspace(0.3,1,1000);
a = 0.5+0.25*N-0.1875*N^2;
b = 1.75-N+0.375*N^2-0.5*(N^2)*CR10.^2;
c = 0.5-N-0.25*N^2+(N^2)*CR10.^2;
p10 = (-b+sqrt(b.^2-4*a*c))/(2*a);
CI10=(N-1)./(2*N*CR10)-p10./(4*CR10);
beta = N/6;
gamma = N^2/48-3*N/64;
delta = N^3/480-3*N^2/256+138*N/3072;
delta0 = beta^2-3*alpha*gamma;
delta1 = -2*beta^3+9*alpha*beta*gamma-27*delta*alpha.^2;
%p9pos = (0.5*(delta1+sqrt(delta1.^2-4*delta0.^3))).^(1/3);
%p9neg = (0.5*(delta1-sqrt(delta1.^2-4*delta0.^3))).^(1/3);
p9pos = 0.5*((N/6)+sqrt((N/6)^2-4*(1-CR9.^2)*((N^2/48)-(3*N/64))))./(1-CR9.^2);
p9neg = 0.5*((N/6)-sqrt((N/6)^2-4*(1-CR9.^2)*((N^2/48)-(3*N/64))))./(1-CR9.^2);
CI9pos = 0.5*N*CR9.*sqrt(0.5*pi*p9pos).*exp(-2*p9pos).*(1-1./(8*p9pos)-9./(128*p9pos.^2)-225./(3072*p9pos.^3));
CI9neg = 0.5*N*CR9.*sqrt(0.5*pi*p9neg).*exp(-2*p9neg).*(1-1./(8*p9neg)-9./(128*p9neg.^2)-225./(3072*p9neg.^3));


x = linspace(0.36,1,1000);
p2 = 0.00512821*(272-sqrt(-66221+505440*x.^2));
M2 = 36*(1-0.5*p2);
z2 = 0.5-0.25*p2;
z2sq = z2.^2; 
gamma2=(1-2*z2sq).*(M2-1)./M2+z2sq.*(2+M2.^2)./(2*M2);
y2 = z2.*(M2-1)./sqrt(2*gamma2.*M2);

%upto p^2,p<<1,b=exp(-p)/2*pi*I0(p),upto (1/N)
N = lambda*timesteps;
xp1 = linspace(0.36,1,1000);
yp1 = 1+(13*xp1-8*xp1.^2-6)./(2*N*xp1.^2);

%upto p^2,p<<1,b=exp(-p)/2*pi*I0(p),upto (1/N^2)
N = lambda*timesteps;
xp2 = linspace(0.36,1,1000);
yp2 = 1+(13*xp2-8*xp2.^2-6)./(72*xp2.^2)-(288*xp2.^3-144*xp2.^4-167*xp2.^2+18)./(4*N^2*xp2.^4);

%upto p^2,p<<1,b=exp(-p)/2*pi*I0(p),upto (1/N^4)
N = lambda*timesteps;
xp4 = linspace(0.36,1,1000);
bracket0 = 1; %Order N^0
bracket1 = (13*xp4-8*xp4.^2-6)./(72*xp4.^2); %Order 1/N
bracket2 = -(288*xp4.^3-144*xp4.^4-167*xp4.^2+18)./(4*N^2*xp4.^4); %Order 1/N^2
bracket3 = (2976*xp4.^5-1216*xp4.^6-2140*xp4.^4+429*xp4.^2-54)./(4*N^3*xp4.^6); %Order 1/N^3
bracket4 = -(245760*xp4.^7-88832*xp4.^8-204640*xp4.^6+60337*xp4.^4-14580*xp4.^2+1620)./(32*N^4*xp4.^8); %Order 1/N^4
yp4 = bracket0+bracket1+bracket2+bracket3+bracket4;


N = lambda*timesteps;
x5 = linspace(0.4,1/sqrt(2),1000);
p5 = 0.5*N*x5.^2./(1+sqrt(1-2*x5.^2));
M=N./(2*p5);
z = sqrt(2*pi)*exp(-2*p5).*(p5.^1.5);
y5 = z.*sqrt(0.5*(M-1));


%using N>>1 approx for line b=f(p),p<<1
N = lambda*timesteps;
x6 = linspace(0.35,1/sqrt(2),1000);

y6 = 1-1./(2*N*x6);

z1(1)=0.6;
z1(2)=0.6;
z2(1)=0.6;
z2(2)=0.5;
z2(3)=0.4;
z2(4)=0.3;




hold off
figure(1); clf
i1 = max(kappa);
i1s = num2str(i1);
i2 = 0.75*max(kappa);
i2s = num2str(i2);
i3 = 0.5*max(kappa);
i3s = num2str(i3);
i4 = 0.25*max(kappa);
i4s = num2str(i4);

sizel(1) = i1;
sizel(2) = i2;
sizel(3) = i3;
sizel(4) = i4;

%for i = 1,length(kappa)
%plot(CR(i),CI(i),'o','MarkerSize',exp(kappa(i)), ...
%     'MarkerFaceColor',kappa1(i),'MarkerEdgeColor',kappa1(i))
%disp(i);
%hold on;
%end



hold off;
figure(1); clf;
fig = scatter(CR,CI,2*exp(kappa),kappa1,'filled'); hold on;
h1 = plot(CR1, CI1, 'k', 'LineWidth', 2); hold on;
h3 = plot(CR3, CI3, 'k', 'LineWidth', 2); hold on;
h4 = plot(CR4, CI4, 'k', 'LineWidth', 2); hold on;
%h13 = plot(x5, y5, 'k--', 'LineWidth', 2); hold on;
h5 = errorbar(CR_exp1,CI_exp1,-dCI_exp1,dCI_exp1,-dCR_exp1,dCR_exp1,'rs','Markerfacecolor','red','MarkerSize',3); hold on;
h6 = errorbar(CR_exp2,CI_exp2,-dCI_exp2,dCI_exp2,-dCR_exp2,dCR_exp2,'rs','Markerfacecolor','red','MarkerSize',6); hold on;
h7 = errorbar(CR_exp3,CI_exp3,-dCI_exp3,dCI_exp3,-dCR_exp3,dCR_exp3,'rs','Markerfacecolor','red','MarkerSize',9); hold on;
h8 = errorbar(CR_exp4,CI_exp4,-dCI_exp4,dCI_exp4,-dCR_exp4,dCR_exp4,'rs','Markerfacecolor','red','MarkerSize',12); hold on;
%h13 = errorbar(CRbar2,CIbar2,-dCR2,dCR2,-dCI2,dCI2,'o','Color',[1 0.5 0.5],'Markeredgecolor',[1 0.5 0.5],'Markerfacecolor',[1 0.5 0.5],'MarkerSize',11);
h13 = scatter(CRbar2,CIbar2,50,'s','cyan','filled'); hold on;
h9 = plot(-z1(1),z2(1),'ks','MarkerSize',10,'Markerfacecolor','black'); hold on;
h10 = plot(-z1(1),z2(2),'ks','MarkerSize',20,'Markerfacecolor','black'); hold on;
h11 = plot(-z1(1),z2(3),'ks','MarkerSize',30,'Markerfacecolor','black'); hold on;
h12 = plot(-z1(1),z2(4),'ks','MarkerSize',40,'Markerfacecolor','black'); 




axis([0 1 -0.3 1]);

s1 = 'b=exp(-p)/2\piI_0(p)';
s2 = 'p=0 parametric';
s3 = 'b=0';
s4 = 'p=0';
s5 = '0 nM/mm';
s6 = '1 nM/mm';
s7 = '5 nM/mm';
s8 = '50 nM/mm';
s9 = strcat('p =',i1s);
s10 = strcat('p =',i2s);
s11 = strcat('p =',i3s);
s12 = strcat('p =',i4s);
s13 = 'approximate analytic line';
s14 = 'p<<1,exp(-M)<<1,order(p^2),order(1/N)';
s15 = 'p<<1,exp(-M)<<1,order(p^2),order(1/N^2)';
s16 = 'p<<1,exp(-M)<<1,order(p^2),order(1/N^4)';


%legend([h5 fig h1 h13],'experiment','simulation','analytic','approximation','Location','northeast');
legend([fig h1 h7 h13],'BPRW theory','BPRW approximation','Experiment','CPM simulation','Location','northeast');



box on;
colormap(jet);

%cb = colorbar; 
%colormap(jet);
%set(cb,'position',[0.21 .175 .045 .4])
xlabel('Chemotactic Ratio');
ylabel('Chemotactic Index');
%set(get(colorbar,'XLabel'),'String','b');
set(gca,'fontsize',15);
saveas(fig,'fig_6_bprw.png');


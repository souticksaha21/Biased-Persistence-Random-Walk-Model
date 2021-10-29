clear all
chifactor = pi;
Z = 1e3; %Number of trial per succesful parameter values
thetab = 0.0*pi; %direction of bias
timesteps = 36; %the total number of timesteps taken
avg_rtime = 1; %the average run time per time step

%alpha/(1-alpha) is the weight given to first/second of the von mises distribution
alphamin = -50; %minimum value of alpha, for this code alpha = 1
alphamax = 50; %maximum value of alpha, for this code alpha = 1
%alpha_pts = 20;

%kappa is the inverse width of one of the von mises distribution
kappamin = 0; %minimum value of kappa
kappamax = 100; %maximum value of kappa
%kappa_pts = 20;

%kappab is the inverse width of the other of the von mises distribution
kappabmin = 0; %minimum value of kappab
kappabmax = 10; %maximum value of kappab
%kappab_pts = 20;

%kappa1 is the weight of the bias distribution
kappa1min = 0; %minimum value of kappa1
kappa1max = 0.2; %maximum value of kappa1
%kappa1_pts = 20;

speed = 10; %speed of the cell

jmax = 100; %no of points for theta sampling in probability
param_pts = 5000; %no of successful parameter points
max_rtime = avg_rtime*timesteps; %total time for which each cell moves

%opening files for writing data
fileID = fopen('bprw_CI_CR_julien.txt','w');
fileID2 = fopen('bprw_full_CI_CR_julien.txt','w');


angle = 0;

%im1 is the label for the total number of successful points sampled
im1 = 0;
%for ia = 1:alpha_pts
%   alpha = alphamin+((alphamax-alphamin)/(alpha_pts-1))*(ia-1);
%  disp(alpha);
% for ik = 1:kappa_pts
%    kappa = kappamin+((kappamax-kappamin)/(kappa_pts-1))*(ik-1);

%   for ikb = 1:kappab_pts
%  kappab = kappabmin+((kappabmax-kappabmin)/(kappab_pts-1))*(ikb-1);
% for ik1 = 1:kappa1_pts
%kappa1 = kappa1min+((kappa1max-kappa1min)/(kappa1_pts-1))*(ik1-1);
while im1 < param_pts
    %alpha = alphamin + (alphamax-alphamin)*rand;
    alpha = 1;
    kappa = kappamin + (kappamax-kappamin)*rand;
    kappa1 = kappa1min + (kappa1max-kappa1min)*rand;
    kappab = kappabmin + (kappabmax-kappabmin)*rand;
    %disp(im1)
    
    %test is the marker for the positivity condition of the probability
    %distribution
    test = 1;
    for j = 1:jmax
        theta = 2*pi*(j-1)/(jmax-1) - pi;
        dtheta = 2*pi/jmax;
        for j1 = 1:jmax
            phi = 2*pi*(j1-1)/(jmax-1) - pi;
            dphi = 2*pi/jmax;
            prob1 = kappa1*cos(theta-thetab);
            prob2 = alpha*exp(kappa*cos(phi))/(2*pi*besseli(0,kappa));
            prob3 = (1-alpha)*exp(kappab*cos(phi))/(2*pi*besseli(0,kappab));
            prob = prob1 + prob2 + prob3;
            %xi=xi+kappa1*cos(theta)*dtheta;
            if prob < 0
                test = 0;
                break
            end
            %psi=psi+
        end
        if test == 0
            break
        end
    end
    if test ~= 0
        disp(im1);
        im1 = im1+1;
        alphach(im1) = alpha;
        kappach(im1) = kappa;
        kappa1ch(im1) = kappa1;
        kappabch(im1) = kappab;
        
    end
end
%end
%end
%end
%end

%disp(im1)
for ik4 = 1:im1
    %initialising chi and psi for analytic calculation
    chi=0;
    psi=0;
    alpha = alphach(ik4);
    kappa = kappach(ik4);
    kappa1 = kappa1ch(ik4);
    kappab = kappabch(ik4);
    disp(ik4)
    for j = 1:jmax
        theta = 2*pi*(j-1)/(jmax-1) - pi; %theta
        dtheta = 2*pi/jmax;
        prob1 = kappa1*cos(theta); %bias distribution
        prob2 = alpha*exp(kappa*cos(theta))/(2*pi*besseli(0,kappa)); %persistence one distribution
        prob3 = (1-alpha)*exp(kappab*cos(theta))/(2*pi*besseli(0,kappab)); %persistence two distribution
        %persistence = prob1+prob2; MISTAKE!!!
        persistence = prob2+prob3;
        bias= prob1;
        chi=chi+bias*cos(theta)*dtheta;
        psi=psi+persistence*cos(theta)*dtheta;
    end
    %chi = chifactor*kappa1;
    %psi = besseli(1,kappa)/besseli(0,kappa);
    lambda=1/avg_rtime;
    lambda_0=lambda*(1-psi);
    bracket1=max_rtime*(1-2*chi^2/(1-psi)^2);
    bracket2=-max_rtime*(chi^2)*exp(-lambda_0*max_rtime)/(1-psi)^2;
    bracket3=((2*chi^2/(1-psi)^2)-1)*(1-exp(-lambda_0*max_rtime))/lambda_0;
    bracket4=lambda_0*(chi^2)*(max_rtime^2)/(2*(1-psi)^2);
    bracket5=(chi^2/(1-psi)^2)*(1-exp(-lambda_0*max_rtime))^2/lambda_0;
    bracket=bracket1+bracket2+bracket3+bracket4+bracket5;
    %analytic CI CR
    CI_analytic=chi*((0.5*lambda/(1-psi))^0.5)*(max_rtime-(1-exp(-lambda_0*max_rtime))/lambda_0)/bracket^0.5;
    CR_analytic=((2/lambda_0)^0.5)*bracket^0.5/max_rtime;
    CR = 0;
    CI = 0;
    distavg = 0;
    distsqravg = 0;
    xavg = 0;
    xsqavg = 0;
    xrratio = 0;
    yavg = 0;
    for i = 1:Z %trial index for each parameter point
        xbac(i) = 0; %initial x position of the cell
        ybac(i) = 0; %initial y position of the cell
        dist = 0; %initial distance the cell covered
        xinit = xbac(i); %initialising cell x position
        yinit = ybac(i); %initialising cell y position
        xold = xbac(i); %storing x position of the cell
        yold = ybac(i); %storing y position of the cell
        thetabac(i) = 2*pi*rand-pi; %randomly initialising velocity direction of the cell
        theta = thetabac(i);
        velxbac(i) = speed*cos(theta); %initial x component of velocity of the cell
        velybac(i) = speed*sin(theta); %initial y component of the velocity of the cell
        for j=1:timesteps
            p=0;
            x=1;
            %choosing run time
            while p < x
                x=rand;
                dtime=max_rtime*rand;
                p=exp(-dtime/avg_rtime)/(avg_rtime*(1-exp(-max_rtime/avg_rtime)));
            end
            dtau=dtime;
            %disp(dtau);
            xbac(i)=xbac(i)+velxbac(i)*dtau; %updating x position
            ybac(i)=ybac(i)+velybac(i)*dtau; %updating y position
            dist=dist+sqrt((xbac(i)-xold)*(xbac(i)-xold)+(ybac(i)-yold)*(ybac(i)-yold)); %updating distance covered
            xold=xbac(i); %saving current x coordinate
            yold=ybac(i); %saving current y coordinate
            thetaold=theta; %saving current theta
            %choosing new theta
            p1=0;
            x1=1;
            while p1 < x1
                x1 = rand;
                thetanew=2*pi*rand-pi;
                phi = thetanew-thetaold;
                prob1 = kappa1*cos(thetanew-thetab);
                prob2 = alpha*exp(kappa*cos(phi))/(2*pi*besseli(0,kappa));
                prob3 = (1-alpha)*exp(kappab*cos(phi))/(2*pi*besseli(0,kappab));
                p1 = prob1 + prob2 + prob3;
                theta=thetanew;
            end
            velxbac(i)=speed*cos(theta);
            velybac(i)=speed*sin(theta);
        end
        distavg = distavg+sqrt((xbac(i)-xinit)*(xbac(i)-xinit)+(ybac(i)-yinit)*(ybac(i)-yinit));
        distsqravg = distsqravg+(xbac(i)-xinit)*(xbac(i)-xinit)+(ybac(i)-yinit)*(ybac(i)-yinit);
        xavg = xavg+(xbac(i)-xinit);
        yavg = yavg+(ybac(i)-yinit);
        xsqavg = xsqavg+(xbac(i)-xinit)*(xbac(i)-xinit);
        xrratio = xrratio+(xbac(i)-xinit)/sqrt((xbac(i)-xinit)*(xbac(i)-xinit)+(ybac(i)-yinit)*(ybac(i)-yinit));
        CR=CR+sqrt((xbac(i)-xinit)*(xbac(i)-xinit)+(ybac(i)-yinit)*(ybac(i)-yinit))/dist;
        deltatot = sqrt((xbac(i)-xinit)*(xbac(i)-xinit)+(ybac(i)-yinit)*(ybac(i)-yinit));
        deltax = xbac(i)-xinit;
        deltay = ybac(i)-yinit;
        cosangle = deltax/deltatot;
        sinangle = deltay/deltatot;
        deltaCI = cosangle*cos(thetab) + sinangle*sin(thetab);
        CI=CI+deltaCI;
    end
    distavg = distavg/Z;
    distsqravg = distsqravg/Z;
    xavg = xavg/Z;
    yavg = yavg/Z;
    xsqavg = xsqavg/Z;
    xrratio = xrratio/Z;
    CI=CI/Z;
    CR=CR/Z;
    x_analytic = speed*(chi/(1-psi))*(max_rtime-(1-exp(-lambda_0*max_rtime))/lambda_0);
    rsqr_analytic = 2*(speed^2)*bracket/lambda_0;
    fprintf(fileID,'%146f,%4.6f,%4.6f,%4.6f,%4.6f,%4.6f,%4.6f,%4.6f\n',x_analytic,rsqr_analytic,distavg,distsqravg,xavg,xsqavg,xrratio,yavg);
    fprintf(fileID2,'%4.6f,%4.6f,%4.6f,%4.6f,%1.6f,%1.6f,%1.6f,%1.6f\n',alpha,kappa,kappab,kappa1,CR,CI,CR_analytic,CI_analytic);
end

fclose(fileID);
fclose(fileID2);




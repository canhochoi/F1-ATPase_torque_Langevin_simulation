%parameters 
global kappa gammaB kBT dtau k01 a1 k02 a2 k21 a21 k03 a3 theta_state method nSim dur
kappa = 56; %pN.nm/rad^2
%gammaG = 3*1e-4; %pN.nm.s
kBT = 4.14; %pN.nm
ktrap = 240; %pN.nm/rad^2
kappaint = 240;

kappain = 223; %pN.nm/rad

kappaout = 150; %pN.nm/rad 
eta = 1e-9; %pN.s/nm^2
% r = 19; %nm
% a = 20; %nm
% gammaB = 8*pi*eta*a^3 + 6*pi*eta*a*r^2; %pN.nm.s
% gammaB = gammaB*2.5;
% gammaG = gammaB;
% gammaint = 3*1e-4;
% gamma = gammaB;
th1 = 3*pi/180; % the first dwell angle in rads (before Pi release)
th2 = 36*pi/180; % the second dwell angle in rads (before ATP binding)
th21 = 72*pi/180;
%th25 = 76*pi/180; % 2.5th state
th3 = 116*pi/180; % the third dwell angle in rads (catalytic dwell)

%theta_state = [th1, th2, th21, th25, th3];
theta_state = [th1, th2, th21, th3];

init_state = [1;0;0;0];

thetaTC = theta_state(init_state == 1) - 5*pi/180;

state = init_state';

%Pi release
% k01 = 6*1.8e4;

%value from previous 
k01 = 5*1.8e+4; % from Noji '10
a1 = 6.7; % from Noji '10 unit is 1/rad 

%ATP binding 
% k02 = 4*8e3;
k02 = 4*9.1e+3; % assumes 0.1 mM ATP concentration
a2 = 0.045*180/pi; % from the stalling experiments unit is 1/rad 

%metastable (related to ADP release)
% k21 = 1/(8e-8);
% k21 = 1/(8e-6);
k21 = 1/(14e-6);  %1/(15e-6) 1/(8e-5) 1/(5e-4) produce 3 lines with dtau = 1e-5
a21 = 0.045*180/pi;

% k25 = 1/(10.e-6);
% a25 = 0.05*180/pi;

%hydrolysis 
a3 = 0;
k03 = 800;

conv = 4184*1e21/(6.02e23); %kcal/mol to pN.nm
% appTorque = -ktrap.*(th-thetaTC);
% dG12 = -2.8*conv-appTorque*(th2-th1); %pN.nm 
% dG23 = -10.2*conv-appTorque*(th3-th2); %pN.nm

%Bronsted slope 
alpha1 = 0.7;
alpha2 = 0.5;
alpha21 = 0.5;
alpha3 = 0;

Pi = 0.01; %mM 
k_01 = Pi*k01/1300;
a_1 = kappa/kBT*(1-alpha1)*(th2-th1);

k_02 = k02/1500;
a_2 = kappa/kBT*(1-alpha2)*(th21-th2);

ADP = 0.01; %mM 
k_21 = ADP*k21/1600;
a_21 = kappa/kBT*(1-alpha21)*(th3-th21);

k_03 = k03/1600;
a_3 = kappa/kBT*(1-alpha3)*(th1+2*pi/3-th3);


%make backward rate large
% k_01 = Pi*k01/1300;
% k_02 = k02*100;
% k_21 = ADP*k21*200;
% k_03 = k03/1600;


% k_01 = 10;
% k_02 = 10;
% k_21 = ADP*k21*200;
% k_03 = k03/1600;

par.k1 =  k01*exp(a1*(-th1)); 
par.k2 =  k02*exp(a2*(-th2));
par.k3 =  k21*exp(a21*(-th21));
par.k4 =  k03*exp(a3*(-th3));

par.k_1 = k_01*exp(-a_1*(-th1)); 
par.k_2 = k_02*exp(-a_2*(-th2)); 
par.k_3 = k_21*exp(-a_21*(-th21)); 
par.k_4 = k_03*exp(-a_3*(-th3)); 

%probability at each state 
prop = @(x,state,par)([par.k1*exp(a1*x).*state(:,1),...  %R1  1 -> 2
                                       par.k2*exp(a2*x).*state(:,2),...  %R2  2 -> 3
                                       par.k3*exp(a21*x).*state(:,3),... %R3  3 -> 4
                                       par.k4*exp(a3*x).*state(:,4),...  %R4  4 -> 1    
                                       par.k_1*exp(-a_1*x).*state(:,2),...  %R5 2 -> 1
                                       par.k_2*exp(-a_2*x).*state(:,3),...  %R6 3 -> 2
                                       par.k_3*exp(-a_21*x).*state(:,4),... %R7 4 -> 3
                                       par.k_4*exp(-a_3*x).*state(:,1),...  %R8 1 -> 4
                                       ]);  
%stochimetry for updating chemical states
%         %R1  R2  R3  R4   R5   R6   R7   R8  
stoch = [-1    0    0   1   1    0   0    -1;... %th1
         1    -1    0   0   -1   1    0     0;...    %th2
         0     1   -1   0   0   -1    1     0;... %th21
         0     0    1  -1   0    0    -1    1;... %th3        
         ]; 

%simulate many trajectories with same duration 
nSim = 50000;

%initial conditions
pos = sqrt(kBT/(10*kappa))*randn(1,nSim)+th1;
% pos_store = NaN*ones(dur,nSim);
% pos_store(1,:) = pos';         
% state_store = NaN*ones(dur,nSim);
% theta_store = NaN*ones(dur,nSim);

%% LG
dur = 120;
dtau = 1e-6; %s
time = (1 : dur)*dtau;

r = 19; %nm
a = 20; %nm
gammaB = 8*pi*eta*a^3 + 6*pi*eta*a*r^2; %pN.nm.s
% gammaB = gammaB*2.5;
gammaB = gammaB*3.36*0.8;

method = 3; %exact based on conditional probability
% [pos_store_m3, state_store_m3, theta_store_m3] = LG(pos,init_state,prop,par,stoch,theta_state);
% [profile_m3, jump_pdf_m3, profile_coarse_m3, jump_pdfcoarse_m3, d_ang_raw_m3, d_ang_rawcoarse_m3, profile_10us_3, jump_10us_3] = angvelplot(pos_store_m3, state_store_m3);

[pos_store_m3_2, state_store_m3_2, theta_store_m3_2] = LG_2(pos,init_state,prop,par,stoch,theta_state);
[profile_m3_2, jump_pdf_m3_2, profile_coarse_m3_2, jump_pdfcoarse_m3_2] = analysis_2(pos_store_m3_2, state_store_m3_2);

%% 3 states 
% make pseudo state very fast (quick check)
% par.k3 =  1e8*exp(a21*(-th21));

%rewrite propensity 
%probability at each state 
prop_3state = @(x,state,par)([par.k1*exp(a1*x).*state(:,1),...  %R1  1 -> 2
                                       par.k2*exp(a2*x).*state(:,2),...  %R2  2 -> 3
                                       par.k4*exp(a3*x).*state(:,3),...  %R4  4 -> 1    
                                       par.k_1*exp(-a_1*x).*state(:,2),...  %R5 2 -> 1
                                       par.k_2*exp(-a_2*x).*state(:,3),...  %R6 3 -> 2
                                       par.k_4*exp(-a_3*x).*state(:,1),...  %R8 1 -> 4
                                       ]);  
%stochimetry for updating chemical states
%         %R1  R2   R4   R5   R6  R8  
stoch_3state = [-1    0    1   1    0   -1;... %th1
         1    -1    0   -1   1    0;...    %th2
         0     1   -1   0    -1    1;... %th3        
         ];             
     
theta_3state = [th1, th2, th3];

init_3state = [1;0;0];

thetaTC_3state = theta_3state(init_3state == 1) - 5*pi/180;

state3 = init_3state';

dur = 120;
dtau = 1e-6; %s
time = (1 : dur)*dtau;

r = 19; %nm
a = 20; %nm
gammaB = 8*pi*eta*a^3 + 6*pi*eta*a*r^2; %pN.nm.s
% gammaB = gammaB*2.5;
gammaB = gammaB*3.36*0.8;

method = 3; %exact based on conditional probability
% [pos_store_m3, state_store_m3, theta_store_m3] = LG(pos,init_state,prop,par,stoch,theta_state);
% [profile_m3, jump_pdf_m3, profile_coarse_m3, jump_pdfcoarse_m3, d_ang_raw_m3, d_ang_rawcoarse_m3, profile_10us_3, jump_10us_3] = angvelplot(pos_store_m3, state_store_m3);

[pos_store_m3_3state, state_store_m3_3state, theta_store_m3_3state] = LG_2(pos,init_3state,prop_3state,par,stoch_3state,theta_3state);
[profile_m3_3state, jump_pdf_m3_3state, profile_coarse_m3_3state, jump_pdfcoarse_m3_3state] = analysis_2(pos_store_m3_3state, state_store_m3_3state);

%%
fig = openfig('averageof4_CI.fig');
axObjs = fig.Children;
dataObjs = axObjs.Children;
ydata = dataObjs(1).YData;
xdata = dataObjs(1).XData;

figure
plot(xdata,ydata,'.','MarkerSize',24)
hold on
plot(profile_coarse_m3_2(:,1),profile_coarse_m3_2(:,2),'LineWidth',2)
plot(profile_coarse_m3_3state(:,1),profile_coarse_m3_3state(:,2),'k--','LineWidth',2)

xlim([-1 121])
legend('Experimental data','4 states','3 states')
xlabel('Angle (deg)')
ylabel('Mean angular jump (deg/ms)')

figure
plot(profile_coarse_m3_2(:,1),profile_coarse_m3_2(:,2))
hold on
plot(profile_m3_2(:,1),profile_m3_2(:,2),'--')
plot(profile_m3_3state(:,1),profile_m3_3state(:,2),'.-')
plot(profile_coarse_m3_3state(:,1),profile_coarse_m3_3state(:,2))

legend('\Deltat=10\mus','\Deltat=1\mus')
xlabel('Angle (deg)')
ylabel('Mean angular jump (deg/ms)')
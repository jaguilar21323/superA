%% X Z M vs u w q dotw
% FLIGHT CONDITION: M = 2 (cruise), h = 45000 ft, rho = 0.0149 lb/ft^3, 
% n = 1, a = 967.7 ft/s

%% Constants
g = ; % gravity at 45000 ft
theta0 = 0; % Euler angle for rotation to final elevation
rho = 0.0149;
cbar = 165; % ft, mean aero dynamic chord
cbarh = 55; % ft
S = 1281; % ft^2, surface area of wing 
Sh = 124.4; % ft^2, surface area of horizontal stab
A = 3.5; % aspect ratio
xw = ; % distance from
e = ; % Oswald's Efficiency Factor
VH = ; % Horizontal Tail Volume Coefficient
swc4 = ; % quarter chord sweep of wing in degrees
a = 967.7; % speed of sound at flight condition
M = 2; % mach at flight condition
W = 80000; % kg, guess for now, have to find what the weight would be around the time cruise starts
u0 = M*a; % velocity at flight condition
q = 0.5*rho*u0^2;
n = 1;
dw = ; % downwash
Cl = n*W/(q*S); % lift coefficient at flight condition
Cd = ; % drag coefficient at flight condition

%% Graph CdvM and xacvM
Cd85 = 0.0105 + 0.2457*Cl^2;
Cd9 = 0.0113 + 0.2387*Cl^2;
Cd95 = 0.0175 + 0.2299*Cl^2;
Cd1p2 = 0.0221 + 0.2128*Cl^2;
Cd1p4 = 0.0243 + 0.2392*Cl^2;
Cd1p6 = 0.0271 + 0.2907*Cl^2;
Cd1p8 = 0.0241 + 0.3636*Cl^2;
Cd2 = 0.0223 + 0.4149*Cl^2;
Cd2p2 = 0.0211 + 0.4717*Cl^2;
Cd = [ Cd85 Cd9 Cd95 Cd1p2 Cd1p4 Cd1p6 Cd1p8 Cd2 Cd2p2 ];
Ma = [ 0.85 0.90 0.95 1.2 1.4 1.6 1.8 2.0 2.2 ];
plot(Ma,Cd)
p = polyfit(M,Cd,3);
k = polyder(p);
CdvM = polyval(k,M);

%% Other Eqts that are for other stuff that we r gonna usel8r oh god plz stop the pain i hate stability 
B = sqrt(1-M^2*cosd(swc4)^2);
Clu = M^2*cosd(swc4)^2*Cl/(1-M^2*cosd(swc4)^2);
ClqW0 = ClalphaW*(0.5+2*xw/cbar);
ClqW = ClqW0*(A+2cosd(swc4))/(A*B+2*cosd(swc4));
ClqH = 2*ClalphaH*nH*VH;
Clq = ClqW+ClqH;
Clalpha = ; % Bens studies (lift curve slope of whole plane)
ClalphaW = ; % lift curve slope of wing
ClalphaH = ; % lift curve slope of horizontal tail
Clalphadot = 2*ClalphaH*nH*VH*dw; % at some point, roskam says to go to Ref. 9 for supersonic methods since this is just for subsonic 
Cwo = (m*g)/(.5*rho*u0^2*S); 
Cdu = Ma * CdvM; % pamadi eqn 4.480 (CdvM = slope of Cd v M)
Cdq = 0;
Cdalphadot = 0;

%% M derivatives
Cmu = -Cl*xacvM*M; % roskam vi eqn 10.12 (xacvM = slope of xac v M)
Cmalpha = -1.4; % from PDR.10
a1 = A*(2(xw/cbar)^2+.5*(xw/cbar))/(A+2*cosd(swc4)); % split Cmqw0 into multiple components bc long
a2 = A^3(tand(swc4)^2)/(24*(A+6*cosd(swc4));
Cmqw0 = -Kw*ClalphaW*cosd(swc4)*(a1+a2+1/8);
Cmwh = -2*(ClalphaH)*nH*VH*(xbarACH-xbarCG);
Cmqw = Cmqw0*((A^3*tand(swc4)^2/(A*B+6*cosd(swc4))+3/B)/(A^3*tand(swc4)^2/(A + 6*cosd(swc4))+3));
Cmq = Cmqw + Cmqh;
Cmalphadot = -2*ClalphaH*nH*VH*(xbarACH-xbarCG)*dw; 

Mu = .5*rho*u0*cbar*Cmu; 
Mw = .5*rho*u0*cbar*Cmalpa; 
Mq = .25*rho*u0*cbar^(-2)*Cmq; 
Mdotw = .25*rho*u0*cbar^(-2)*Cmalphadot; 


%% Z derivatives this is where i like to edit things 
Czu = -2 * Cl - Clu ; % pamadi eqn 4.439
Czalpha = -Clalpha; 
Czq = -Clq;
Czalphadot = -Clalphadot; % pamadi eqn 4.457

Zu = -rho*u0*S*Cw0*sin(theta0) + .5*rho*u0*S*Czu;
Zw = .5*rho*u0*S*Czalpha; 
Zq = .25*rho*u0*cbar*S*Czq; 
Zwdot = .25*rho*cbar*S*Czalphadot; 

%% X derivatives 
Cxu = -2 * Cd - Cdu; % eqn 4.432 pamadi
Cxalpha = Cl-2*Cl*Clalpha/(pi*S*e); % etkin eqt 5.2,2
Cxq = -Cdq; % pamadi eqn 4.456
Cxdotalpha = -Cdalphadot; % pamadi eqn 4.455

Xu = rho*u0*S*Cw0*sin(theta0) + .5*rho*u0*S*Cxu; 
Xw = .5*rho*u0*S*Cxalpha; 
Xq = .25*rho*u0*cbar*S*Cxq; 
Xdotw = .25*rho*cbar*S*Cxalphadot; 



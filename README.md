# Brouwer
Brouwer's Theory on Artificial Satellite Code
close all
clear all
clc

RE = 6378.137e3;                                                                % mean radius of Earth. meters
mu = 398600.4418e9;                                                         % gravitational parameter. m^3/s^2
%%%%%%%% UT at which the orbital parameters are recorded.
ee = .0977;                                                               % eccentricity at epoch of observation
i = deg2rad(51.6203);                                                       % inclination at epoch of observation. radians
% r_p = 397e3+RE;                                                             % radius of perigee. meters
% r_a = 411e3+RE;                                                             % radius of apogee. meters
r_p = 35786e3+RE;                                                             % radius of perigee. meters
r_a = 35800e3+RE;                                                             % radius of apogee. meters
Omega = deg2rad(353.5939);                                                  % right ascension of ascending node.
omega = deg2rad(198.2108);                                                  % argument of perigee.
aa = (r_p+r_a)/2;                                                            % semi major axis of the orbit. meters
Period = 2*pi*sqrt(aa^3/mu);                                                 % orbit period. seconds
p = aa*(1-ee^2);                                                              % semi-latus rectum. meters

thetastar = linspace(0,2*pi,1e3);                                           % varying true anomaly from 0deg to 360 deg for one complete orbit data
r = p./(1+ee.*cos(thetastar));                                               % magnitude of distance from Earth center to ISS at each value of true anomaly
theta = thetastar + omega;                                                  % true anomaly as measured from line of nodes. radians
DCM = zeros(3,3,length(thetastar));                                         % memory allocation for direction cosine matrix
pos_rTh = zeros(1,3,length(thetastar));                                     % memory allocation for storing ISS position ('r-theta-h' frame). meters
pos_ECI = zeros(1,3,length(thetastar));                                     % memory allocation for storing ISS position (ECI frame). meters
for j = 1:length(thetastar)
    pos_rTH(:,:,j)=r(j)*[1 0 0];                                            % position vector of ISS in 'r-theta-h' frame. meters
    DCM(:,:,j) = [cos(Omega)*cos(theta(j))-sin(Omega)*cos(i)*sin(theta(j)) -cos(Omega)*sin(theta(j))-sin(Omega)*cos(i)*cos(theta(j)) sin(Omega)*sin(i);...
        sin(Omega)*cos(theta(j))+cos(Omega)*cos(i)*sin(theta(j)) -sin(Omega)*sin(theta(j))+cos(Omega)*cos(i)*cos(theta(j)) -cos(Omega)*sin(i); ...
        sin(i)*sin(theta(j)) sin(i)*cos(theta(j)) cos(i)];
    pos_ECI(:,:,j) = pos_rTH(:,:,j)*DCM(:,:,j)';                            % position vector of ISS in ECI frame. meters
end
figure(1)
plot3(squeeze(pos_ECI(1,1,:)./1e3),squeeze(pos_ECI(1,2,:)./1e3),squeeze(pos_ECI(1,3,:)./1e3))
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
title('ECI Orbits')
grid on
hold on


e_dd = ee;
f_d = 0;
r_d = r_p;
I_dd = i;
J2 = 1.0826e-3;
J3 = -2.532e-6;
J4 = -1.649e-6;
J5 = -.210e-6;
k2 = J2*RE^2/2;
A_30 = -J3*RE^3;
k4 = -(3*J4/8)*RE^4;
A_50 = -J5*RE^5;
a_dd = aa;
t = 0;
g0_dd = omega;
h0_dd = Omega;


counter = 1;
while (1==1)
    eta = sqrt(1-e_dd^2);
    theta = cos(I_dd);
    gamma_2 = k2/a_dd^2;
    gamma_4 = k4/a_dd^4;
    gamma_3 = A_30/a_dd^3;
    gamma_5 = A_50/a_dd^5;
    gamma_2_d = gamma_2*(eta^(-4));
    gamma_4_d = gamma_4*(eta^(-8));
    gamma_3_d = gamma_3*(eta^(-6));
    gamma_5_d = gamma_5*(eta^(-10));
    
    l0_dd = sqrt(mu/a_dd^3)*t;
    
    n0 = sqrt(mu/a_dd^3);
    l_dd = n0*t*(1+1.5*gamma_2_d*eta*(-1+3*theta^2)+(3/32)*gamma_2_d^2*eta*(-15+16*eta+25*eta^2+(30-96*eta-90*eta^2)*theta^2 +(105+144*eta+25*eta^2)*theta^4 ...
        )+(15/16)*gamma_4_d*eta*e_dd^2*(3-30*theta^2+35*theta^4) ) + l0_dd;
    g_dd = n0*t*(1.5*gamma_2_d*(-1+5*theta^2)+(3/32)*gamma_2_d^2*(-35+24*eta+25*eta^2+(90-192*eta-126*eta^2)*theta^2+(385*+360*eta+45*eta^2)*theta^4)+ ...
        (5/16)*gamma_4_d*(21-9*eta^2+(-270+126*eta^2)*theta^2+(385-189*eta^2)*theta^4)) + g0_dd;
    h_dd = n0*t*(-3*gamma_2_d*theta+(3/8)*gamma_2_d^2*((-5+12*eta+9*eta^2)*theta+(-35-36*eta-5*eta^2)*theta^3)+(5/4)*gamma_4_d*(5-3*eta^2)*theta*(3-7*theta^2) )+h0_dd;
    
    
    delta1_e = ((1/8)*gamma_2_d*e_dd*eta^2*(1-11*theta^2-40*theta^4*(1-5*theta^2)^(-1))-(5/12)*(gamma_4_d/gamma_2_d)*e_dd*eta^2*(1-3*theta^2-8*theta^4*(1-5*theta^2)^(-1)) )*cos(2*g_dd) ...
        + ((1/4)*(gamma_3_d/gamma_2_d)*eta^2*sin(I_dd)+(5/64)*(gamma_5_d/(gamma_2_d*eta^2))*sin(I_dd)*(4+3*e_dd^2)*(1-9*theta^2-24*theta^4*(1-5*theta^2)^(-1)) )*sin(g_dd) ...
        -(35/384)*(gamma_5_d/gamma_2_d)*e_dd^2*eta^2*sin(I_dd)*(1-5*theta^2-16*theta^4*(1-5*theta^2)^(-1))*sin(3*g_dd);
    delta1_I = -e_dd*delta1_e/(eta^2*tan(I_dd));
    l_d = l_dd + ((1/8)*gamma_2_d*eta^3*(1-11*theta^2-40*theta^4*(1-5*theta^2)^(-1)) - (5/12)*(gamma_4_d/gamma_2_d)*eta^3*(1-3*theta^2-8*theta^4*(1-5*theta^2)^(-1)) )*sin(2*g_dd) ...
        + (-(1/4)*(gamma_3_d/gamma_2_d)*(eta^3/e_dd)*sin(I_dd)-(5/64)*(gamma_5_d/gamma_2_d)*(eta^3/e_dd)*sin(I_dd)*(4+9*e_dd^2)*(1-9*theta^2-24*theta^4*(1-5*theta^2)^(-1)) )*cos(g_dd) ...
        + (35/384)*(gamma_5_d/gamma_2_d)*eta^3*e_dd*sin(I_dd)*(1-5*theta^2-16*theta^4*(1-5*theta^2)^(-1))*cos(3*g_dd);
    g_d = g_dd + (-(1/16)*gamma_2_d*((2+e_dd^2)-11*(2+3*e_dd^2)*theta^2-40*(2+5*e_dd^2)*theta^4*(1-5*theta^2)^(-1)-400*e_dd^2*theta^6*(1-5*theta^2)^(-2)) ...
        +(5/24)*(gamma_4_d/gamma_2_d)*(2+e_dd^2-3*(2+3*e_dd^2)*theta^2-8*(2+5*e_dd^2)*theta^4*(1-5*theta^2)^(-1)-80*e_dd^2*theta^6*(1-5*theta^2)^(-2)))*sin(2*g_dd) ...
        + ((1/4)*(gamma_3_d/gamma_2_d)*(sin(I_dd)/e_dd - ((e_dd*theta^2)/(sin(I_dd))))+(5/64)*(gamma_5_d/gamma_2_d)*(((eta^2*sin(I_dd))/e_dd - e_dd*theta^2/ ...
        sin(I_dd))*(4+3*e_dd^2)+e_dd*sin(I_dd)*(26+9*e_dd^2))*(1-9*theta^2-24*theta^4*(1-5*theta^2)^(-1))-(15/32)*(gamma_5_d/gamma_2_d)*e_dd*theta^2*sin(I_dd)* ...
        (4+3*e_dd^2)*(3+16*theta^2*(1-5*theta^2)^(-1)+40*theta^4*(1-5*theta^2)^(-2)) )*cos(g_dd) + ((-35/1152)*(gamma_5_d/gamma_2_d)*(e_dd*sin(I_dd)*(3+2*e_dd^2)- ...
        e_dd^3*theta^2/sin(I_dd))*(1-5*theta^2-16*theta^4*(1-5*theta^2)^(-1))+(35/576)*(gamma_5_d/gamma_2_d)*e_dd^3*theta^2*sin(I_dd)*(5+32*theta^2*(1-5*theta^2)^(-1)+80*theta^4*(1-5*theta^2)^(-2)))*cos(3*g_dd);
    h_d = h_dd + ((-1/8)*gamma_2_d*e_dd^2*theta*(11+80*theta^2*(1-5*theta^2)^(-1)+200*theta^4*(1-5*theta^2)^(-2))+(5/12)*(gamma_4_d/gamma_2_d)*e_dd^2*...
        theta*(3+16*theta^2*(1-5*theta^2)^(-1)+40*theta^4*(1-5*theta^2)^(-2)) )*sin(2*g_dd) + ((1/4)*(gamma_3_d/gamma_2_d)*(e_dd*theta/sin(I_dd))+(5/64)* ...
        (gamma_5_d/gamma_2_d)*(e_dd*theta/sin(I_dd))*(4+3*e_dd^2)*(1-9*theta^2-24*theta^4*(1-5*theta^2)^(-1))+(15/32)*(gamma_5_d/gamma_2_d)*e_dd*theta*sin(I_dd)* ...
        (4+3*e_dd^2)*(3+16*theta^2*(1-5*theta^2)^(-1)+40*theta^4*(1-5*theta^2)^(-2)) )*cos(g_dd) + ((-35/1152)*(gamma_5_d/gamma_2_d)*(e_dd^3*theta/sin(I_dd)) ...
        *(1-5*theta^2-16*theta^4*(1-5*theta^2)^(-1))-(35/576)*(gamma_5_d/gamma_2_d)*e_dd^3*theta*sin(I_dd)*(5+32*theta^2*(1-5*theta^2)^(-1)+80*theta^4*(1-5*theta^2)^(-2)) )*cos(3*g_dd);
    
    a = a_dd*(1+gamma_2*((-1+3*theta^2)*(a_dd^3/r_d^3 - eta^(-3))+3*(1-theta^2)*(a_dd^3/r_d^3)*cos(2*g_d+2*f_d)));
    e = e_dd + delta1_e+ (eta^2/(2*e_dd))*(gamma_2*((-1+3*theta^2)*(a_dd^3/r_d^3 - eta^(-3))+3*(1-theta^2)*(a_dd^3/r_d^3 - eta^(-4))*cos(2*g_d+2*f_d))- ...
    gamma_2_d*(1-theta^2)*(3*e_dd*cos(2*g_d+f_d)+e_dd*cos(2*g_d+3*f_d)) );
    I = I_dd + delta1_I + .5*gamma_2_d*theta*sqrt(1-theta^2)*(3*cos(2*g_d+2*f_d)+3*e_dd*cos(2*g_d+f_d)+e_dd*cos(2*g_d+3*f_d));
    l = l_d - (eta^3/(4*e_dd))*gamma_2_d*(2*(-1+3*theta^2)*(a_dd^2*eta^2/r_d^2 + a_dd/r_d + 1)*sin(f_d) + 3*(1-theta^2)*((-a_dd^2*eta^2/r_d^2 - a_dd/r_d +1)*sin(2*g_d+f_d)+(a_dd^2*eta^2/r_d^2 + a_dd/r_d + 1/3)*sin(2*g_d+3*f_d)) );
    g = g_d + (eta^2/(4*e_dd))*gamma_2_d*(2*(-1+3*theta^2)*(a_dd^2*eta^2/r_d^2 + a_dd/r_d + 1)*sin(f_d) + 3*(1-theta^2)*((-a_dd^2*eta^2/r_d^2 - a_dd/r_d +1)*sin(2*g_d+f_d)+(a_dd^2*eta^2/r_d^2 + a_dd/r_d + 1/3)*sin(2*g_d+3*f_d)) )...
        +(1/4)*gamma_2_d*(6*(-1+5*theta^2)*(f_d-l_d+e_dd*sin(f_d))+(3-5*theta^2)*(3*sin(2*g_d+2*f_d)+3*e_dd*sin(2*g_d+f_d)+e_dd*sin(2*g_d+3*f_d)) );
    h = h_d - .5*gamma_2_d*theta*(6*(f_d-l_d+e_dd*sin(f_d))-3*sin(2*g_d+2*f_d)-3*e_dd*sin(2*g_d+f_d)-e_dd*sin(2*g_d+3*f_d));
    
    options = optimset('Display','off');
    [E_d,expn] = fsolve(@(x) [x-e*sin(x)-l],l,options);
    r_d = a*(1-e*cos(E_d));
    f_d = acos(a*(cos(E_d)-e)/r_d);
    if (l>pi )
        f_d = 2*pi-f_d;
    end
    x = r_d*(cos(g+f_d)*cos(h)-sin(g+f_d)*sin(h)*cos(I));
    y = r_d*(cos(g+f_d)*sin(h)+sin(g+f_d)*cos(h)*cos(I));
    z = r_d*sin(g+f_d)*sin(I);
    
    x_per(counter) = x;
    y_per(counter) = y;
    z_per(counter) = z;
    
    e_dd = e;
    I_dd = I;
    a_dd = a;
    t = t+Period/500;
    g0_dd = g;
    h0_dd = h;
    if (t > Period/2)
    break;
    end
    counter = counter + 1;
end
figure(1)
plot3(x_per./1e3,y_per./1e3,z_per./1e3,'r*-')
legend('No Higher Order Terms','Brouwer')



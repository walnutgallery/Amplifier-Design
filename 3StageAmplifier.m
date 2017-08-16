clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%% 3st STAGE %%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%% 3st STAGE %%%%%%%%%%%%%%%%%%%%%%%%%
syms Rl_3 Cl_3 Vc_3 Rc_3 Re_3 Rb1_3 Rb2_3 Rs_3 Cl_3 Cb vi_3 re_3 beta_3 
syms Icq_3 Vth_3 Vcc alpha_3 vin1 re_3 Rb_3
syms rpi_3 Rs
Rl_3 = 16;
beta_3 = 110;
Vth_3 = 26/1000;
Vcc = 20;
alpha_3 = beta_3/(beta_3+1);


Rb1_3 = 20000;
Rb2_3 = 20000;
Re_3 = 100;
Rb_3 = parallel(Rb1_3, Rb2_3);

Rs = 0;
VbbEquation1 =Rb2_3/(Rb1_3+Rb2_3)*Vcc;
VbbEquation2 = .7 +(Rb_3/beta_3+((1+beta_3)/beta_3)*Re_3)*Icq_3;
solvedIcq = vpa(solve(VbbEquation1 == VbbEquation2, Icq_3));
Ieq_3=solvedIcq/alpha_3;
rpi_3=(beta_3+1)*Vth_3/Ieq_3;

re_3= rpi_3/(beta_3+1);
R3ParR1 = parallel(Re_3,Rl_3);
temp2 =rpi_3+(1+beta_3)*R3ParR1;
zi_3 = parallel(temp2,Rb_3)
%pretty(zi_3)
Av3=zi_3/(Rs+zi_3)*(R3ParR1/(re_3+R3ParR1))
%pretty(Av3)


%%%%%%%%%%%%%%%%%%%%%%%%% 2nd STAGE %%%%%%%%%%%%%%%%%%%%%%%%%
syms Rl_2 Cl_2 Vc_2 Rc_2 Re_2 Rb1_2 Rb2_2 Rs_2 Cl_2 Cb vi_2 re_2 beta_2 
syms Icq_2 Vth_2 Vcc alpha_2 vin1

Vth_2 = 26/1000; %given
beta_2 = 255; %given
alpha_2= beta_2/(beta_2+1); %given
vcpk = Av*vin; %picked based on Av
Vopk = vcpk

vin = 25/1000; %picked
Vcc = 10; %picked
Rb1_2 = 1500000; %picked
Rb2_2 = 500000; %picked
Av = 30; %picked
vbepk = 5/1000; %picked
Zi=10000;
Rl_2 = zi_3; %chosen from input impedance of the common collector
Vcemin_2 = .5;
%wanted
%Av = 30
%zin = small
%large gain
%Need Re and Rc
%Vcc = Vce + Icq*(Rc+(Re1+Re2)/alpha);

%equations to solve for Re and Rc
Re_2 = (Zi/(1+beta_2))-alpha_2*Vth_2/Icq_2; %equation 5.65a
%IcqEquation1_2 = (Vcc-vcpk*(1+1/Av)-.5)/(Rc_2+Zi/beta_2); %equation 5.66
IcqEquation1_2 = Vth/(parallel(Rc_2, Rl_2))*vin/vbepk; %equation 5.72
IcqEquation2_2 = beta_2*Vth_2/(Zi)*(vin/vbepk-1) %equation 5.73
Icq_2 = IcqEquation2_2
solvedRc_2 = vpa(solve(IcqEquation1_2==IcqEquation2_2, Rc_2))
solvedRe_2 = subs((Zi/(1+beta_2))-alpha_2*Vth_2/Icq_2, Icq_2, IcqEquation2_2)
re_2 = vpa(Vth_2/IcqEquation2_2)
%actual calculations
 Rb_2=Rb1_2*Rb2_2/(Rb1_2+Rb2_2)
% 
 zi = (re_2 + solvedRe_2)*(1 + beta_2)*Rb_2/(Rb_2 + (re_2 + solvedRe_2)*(1 + beta_2)*Rb_2)
% actualAv = (Icq_2*Rc_2/Vth_2)/(1+Icq_2*Re_2/Vth_2);
% zo = Rc_2;


%%%%%%%%%%%%%%%%%%%%%%%%% 1st STAGE %%%%%%%%%%%%%%%%%%%%%%%%%


syms Rl_1 Cl_1 Vc_1 Rc_1 Re_1 Rb1_1 Rb2_1 Rs_1 Cl_1 Cb vi_1 re_1 beta_1 
syms Icq_1 Vth_1 Vcc alpha_1 vin1
vin = 10/1000; %picked
Vcc = 10; %picked
Rb1_1 = 200000; %picked
Rb2_1 = 300000; %picked
Vth_1 = 26/1000; %given
beta_1 = 255; %given
alpha_1= beta_1/(beta_1+1); %given
vcpk = 25/1000; %picked based on Av
Av = 2.5; %picked
vbepk = 5/1000; %picked
Zi=200000;
%wanted
%Av = 2-3 for this 2.5
%zin = 200k
%limited gain large impedance
%Need Re and Rc
%Vcc = Vce + Icq*(Rc+(Re1+Re2)/alpha);

%equations to solve for Re and Rc
Re_1 = (Zi/(1+beta_1))-alpha_1*Vth_1/Icq_1; %equation 5.65a
%IcqEquation1_1 = (Vcc-vcpk*(1+1/Av)-.5)/(Rc_1+Zi/beta_1); %equation 5.66
IcqEquation1_1 = Vth/(parallel(Rc_1, Rl_1))*vin/vbepk; %equation 5.72
IcqEquation2_1 = beta_1*Vth_1/(Zi)*(vin/vbepk-1) %equation 5.73
Icq_1 = IcqEquation2_1
solvedRc_1 = vpa(solve(IcqEquation1_1==IcqEquation2_1, Rc_1))
solvedRe_1 = subs((Zi/(1+beta_1))-alpha_1*Vth_1/Icq_1, Icq_1, IcqEquation2_1)
re_1 = vpa(Vth_1/IcqEquation2_1)
%actual calculations
 Rb_1=Rb1_1*Rb2_1/(Rb1_1+Rb2_1)
% 
 zi = (re_1 + solvedRe_1)*(1 + beta_1)*Rb_1/(Rb_1 + (re_1 + solvedRe_1)*(1 + beta_1)*Rb_1)
% actualAv = (Icq_1*Rc_1/Vth_1)/(1+Icq_1*Re_1/Vth_1);
% zo = Rc_1;

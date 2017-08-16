clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%% 3st STAGE %%%%%%%%%%%%%%%%%%%%%%%%%
syms Rl_3 Cl_3 Vc_3 Rc_3 Re_3 Rb1_3 Rb2_3 Rs_3 Cl_3 Cb vi_3 re_3 beta_3 
syms Icq_3 Vth_3 Vcc alpha_3 vin1 re_3 Rb_3
syms rpi_3 Rs
Rl_3 = 16;
beta_3 = 376;
Vth_3 = 26/1000;
Vcc = 10;
alpha_3 = beta_3/(beta_3+1);


Rb1_3 = 500000;
Rb2_3 = 1000000;
Re_3 = 50;
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


syms Rl_2 Cl_2 Vc_2 Rc_2 Re_2 Rb1_2 Rb2_2 Rs_2 Cl_2 Cb vi_2 re_2 beta_2 
syms Icq_2 Vth_2 Vcc alpha_2 vin1


vin = 25/1000; %picked
Vcc = 10; %picked
Rb1_2 = 200000; %picked
Rb2_2 = 400000; %picked
Av_2 = 30; %picked
vbepk = 5/1000; %picked
Zi=10000;
Rl_2 = zi_3; %chosen from input impedance of the common collector
Vcemin_2 = .5;


Vth_2 = 26/1000; %given
beta_2 = 255; %given
alpha_2= beta_2/(beta_2+1); %given
vcpk = Av_2*vin; %picked based on Av
Vopk = vcpk


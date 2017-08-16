clear; clc; close all;
%declare equations
syms beta icq vth rc rl vcc vcpk av rc zi vipk vbepk re Re Rb;
vcc = 10; %assumed
zi = 190000; %picked
vbepk = 5/1000; %assumed
beta = 255; %assumed
av = 2.6; %picked
vth = .026; %constant
vipk = .2; %picked
rl=    9.6154e+03; %picked from 2nd stage
vcpk = vipk*abs(av); %picked

Rb1 = 5000000;
Rb2 = 5000000;
alpha=beta/(beta+1);

%icq equations
icqEquation = (vcc - vcpk*(1+1/av)-.5)/(rc+zi/beta);
icqEquation2 = vth/(rc*rl/(rc+rl))*(abs(vipk/vbepk))*abs(av);
icqEquation3=((beta*vth)/zi)*((abs(vipk/vbepk))-1);

%solve for the resistor values
solvedRc = vpa(solve(icqEquation2==icqEquation3,rc))

solvedRe = zi/(1+beta)-vth/icqEquation3

solvedre = vth/icqEquation3

Rb = parallel(Rb1,Rb2);
Rs = 0;
rpi=(beta+1)*solvedre;
zi=parallel((rpi+(1+beta)*(solvedRe)),Rb)
Rl=1371.889;
%Av = -(zi/(Rs+zi))*(alpha*parallel(solvedRc,Rl))/(solvedre+parallel(solvedRc,Rl))

clear; clc; close all;
%declare equations
syms beta icq vth rc rl vcc vcpk av rc zi vipk vbepk re Re Rb;
vcc = 10; %assumed
zi = 15000; %picked
vbepk = 5/1000; %assumed
beta = 255; %assumed
av = -10; %picked
vth = .026; %constant
vipk = .25; %picked
rl= 1200; %picked
vcpk = vipk*abs(av); %picked

Rb1 = ;
Rb2 = 3000000;
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
Rl=1800;
%Av = -(zi/(Rs+zi))*(alpha*parallel(solvedRc,Rl))/(solvedre+parallel(solvedRc,Rl))

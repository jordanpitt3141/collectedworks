% Plot the results from the second-order SWWE solver
% Segar_SWWE_rk_2_2_SWWE.for Figure 2

clear all;
close all;

B = load('Segur_rk_2_2_SWWE.r');
load D:\f77l\Boussq\Segur\Fig2-joeharvey78\Fig2a.txt
load D:\f77l\Boussq\Segur\Fig2-joeharvey78\Fig2b.txt
load D:\f77l\Boussq\Segur\Fig2-joeharvey78\Fig2c.txt
load D:\f77l\Boussq\Segur\Fig2-joeharvey78\Fig2d.txt
load D:\f77l\Boussq\Segur\Fig2-joeharvey78\Fig2e.txt

tt    = B(:,1);
h000  = B(:,2);
h050  = B(:,3);
h100  = B(:,4);
h150  = B(:,5);
h200  = B(:,6);

g = 9.81;
h0 = 0.1;
L = 31.6;
b = 2*0.61;

figure(1)
plot(tt*sqrt(g/h0),1.5*(h000 - h0)/h0,'or','MarkerSize',2)
%axes('FontSize',16,'FontName','Times','Box','on')
%axis([0 40 0.18 .24])
%hold on
%text(10,.23,'x = 10.8m','FontAngle','italic','FontSize',20,'FontName','Times')
%xlabel('t(s)','FontAngle','italic','FontSize',20,'FontName','Times')
%ylabel('h(m)','FontAngle','italic','FontSize',20,'FontName','Times')
%axis square
title('SWWE Segur x/h0 = 0')
xlabel('t*sqrt{g/h0} - x/h0')
ylabel('(h - h0)/h0')
xlim([-10 250])
ylim([-0.1001 0.1001])
hold on
plot(Fig2a(:,1),Fig2a(:,2))

figure(2)
plot(tt*sqrt(g/h0) - 50,1.5*(h050-h0)/h0,'or','MarkerSize',2)
title('SWWE Segur x/h = 50')
xlabel('t*sqrt{g/h} - x/h')
ylabel('(h - h0)/h0')
xlim([-10 250])
ylim([-0.1001 0.1001])
hold on
plot(Fig2b(:,1),Fig2b(:,2));


figure(3)
plot(tt*sqrt(g/h0) - 100,1.5*(h100-h0)/h0,'or','MarkerSize',2)
title('SWWE Segur x/h = 100')
xlabel('t*sqrt{g/h} - x/h')
ylabel('(h - h0)/h0)')
xlim([-10 250])
ylim([-0.1001 0.1001])
hold on
plot(Fig2c(:,1),Fig2c(:,2));

figure(4)
plot(tt*sqrt(g/h0) - 150,1.5*(h150-h0)/h0,'or','MarkerSize',2)
title('SWWE Segur x/h = 150')
xlabel('t*sqrt{g/h} + x/h')
ylabel('(h - h0)/h0')
xlim([-10 250])
ylim([-0.1001 0.1001])
hold on
plot(Fig2d(:,1),Fig2d(:,2));

figure(5)
plot(tt*sqrt(g/h0) - 200,1.5*(h200-h0)/h0,'or','MarkerSize',2)
title('SWWE Segur x/h = 200')
xlabel('t*sqrt{g/h} - x/h')
ylabel('(h - h0)/h0')
xlim([-10 250])
ylim([-0.1001 0.1001])
hold on
plot(Fig2e(:,1),Fig2e(:,2));

% Initial and final profile
A = load('Segur_rk_2_2_SWWE.out');
n = max(size(A));

x  = A(1:n/2,1);

% water depth
figure(6)
subplot(2,1,1)
plot(x,A(1:n/2,2))
title('t = 0s')
xlabel('x(m)')
ylabel('h(m)')
xlim([0 100])
ylim([0.06 0.11])
subplot(2,1,2)
plot(x,A(n/2+1:n,2))
title('t = 23s')
xlabel('x(m)')
ylabel('h(m)')
xlim([0 100])
ylim([0.06 0.11])

% water velocity
figure(7)
subplot(2,1,1)
plot(x,A(1:n/2,3))
title('t = 0s')
xlabel('x')
ylabel('u(m/s)')
xlim([0 60])
ylim([-0.04 0.04])
subplot(2,1,2)
plot(x,A(n/2+1:n,3))
title('t = 23s')
xlabel('x')
ylabel('u(m/s)')
xlim([0 100])
ylim([-0.08 0.08])

figure(8)
plot(tt,h000,'or','MarkerSize',2)
%axes('FontSize',16,'FontName','Times','Box','on')
%axis([0 40 0.18 .24])
%hold on
%text(10,.23,'x = 10.8m','FontAngle','italic','FontSize',20,'FontName','Times')
%xlabel('t(s)','FontAngle','italic','FontSize',20,'FontName','Times')
%ylabel('h(m)','FontAngle','italic','FontSize',20,'FontName','Times')
%axis square
title('SWWE Segur x = 50.5m')
xlabel('t(s)')
ylabel('h(m)')
xlim([0 50])
ylim([0.08 0.12])

figure(9)
plot(tt,h050,'or','MarkerSize',2)
title('SWWE Segur x = 55m')
xlabel('t(s)')
ylabel('h(m)')
xlim([0 50])
ylim([0.08 0.12])

figure(10)
plot(tt,h100,'or','MarkerSize',2)
title('SWWE Segur x = 60m')
xlabel('t(s)')
ylabel('h(m)')
xlim([0 50])
ylim([0.08 0.12])

figure(11)
plot(tt,h150,'or','MarkerSize',2)
title('SWWE Segur x = 65m')
xlabel('t(s)')
ylabel('h(m)')
xlim([0 50])
ylim([0.08 0.12])

figure(12)
plot(tt,h200,'or','MarkerSize',2)
title('SWWE Segur x = 70m')
xlabel('t(s)')
ylabel('h(m)')
xlim([0 50])
ylim([0.08 0.12])

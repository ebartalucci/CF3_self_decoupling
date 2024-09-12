load CF3_15_90_100_C00.mat
load CF3_15_90_100_C10.mat
load CF3_15_90_100_C20.mat
load CF3_15_90_100_C25.mat
load CF3_15_90_100_C30.mat
load CF3_15_90_100_C40.mat
load CF3_15_90_100_C50.mat
load CF3_15_90_100_C55.mat
load CF3_15_90_100_C60.mat
load CF3_15_90_100_C70.mat
load CF3_15_90_100_C80.mat
load CF3_15_90_100_C90.mat

load CF3_15_90_312_C00.mat
load CF3_15_90_312_C10.mat
load CF3_15_90_312_C20.mat
load CF3_15_90_312_C25.mat
load CF3_15_90_312_C30.mat
load CF3_15_90_312_C40.mat
load CF3_15_90_312_C50.mat
load CF3_15_90_312_C55.mat
load CF3_15_90_312_C60.mat
load CF3_15_90_312_C70.mat
load CF3_15_90_312_C80.mat
load CF3_15_90_312_C90.mat

load CF3_15_90_420_C00.mat
load CF3_15_90_420_C10.mat
load CF3_15_90_420_C20.mat
load CF3_15_90_420_C25.mat
load CF3_15_90_420_C30.mat
load CF3_15_90_420_C40.mat
load CF3_15_90_420_C50.mat
load CF3_15_90_420_C55.mat
load CF3_15_90_420_C60.mat
load CF3_15_90_420_C70.mat
load CF3_15_90_420_C80.mat
load CF3_15_90_420_C90.mat

load CF3_15_90_160_C00.mat
load CF3_15_90_160_C10.mat
load CF3_15_90_160_C20.mat
load CF3_15_90_160_C25.mat
load CF3_15_90_160_C30.mat
load CF3_15_90_160_C40.mat
load CF3_15_90_160_C50.mat
load CF3_15_90_160_C55.mat
load CF3_15_90_160_C60.mat
load CF3_15_90_160_C70.mat
load CF3_15_90_160_C80.mat
load CF3_15_90_160_C90.mat

load CF3_15_90_670_C00.mat
load CF3_15_90_670_C10.mat
load CF3_15_90_670_C20.mat
load CF3_15_90_670_C25.mat
load CF3_15_90_670_C30.mat
load CF3_15_90_670_C40.mat
load CF3_15_90_670_C50.mat
load CF3_15_90_670_C55.mat
load CF3_15_90_670_C60.mat
load CF3_15_90_670_C70.mat
load CF3_15_90_670_C80.mat
load CF3_15_90_670_C90.mat

data100A=zeros(12,32768);
data312A=zeros(12,32768);
data420A=zeros(12,32768);
data160A=zeros(12,32768);
data670A=zeros(12,32768);

data100A(01,:)=phase2(CF3_15_90_100_C00,15000,32768,130,0,0,0,0);
data100A(02,:)=phase2(CF3_15_90_100_C10,15000,32768,130,0,0,0,0);
data100A(03,:)=phase2(CF3_15_90_100_C20,15000,32768,130,0,0,0,0);
data100A(04,:)=phase2(CF3_15_90_100_C25,15000,32768,130,0,0,0,0);
data100A(05,:)=phase2(CF3_15_90_100_C30,15000,32768,130,0,0,0,0);
data100A(06,:)=phase2(CF3_15_90_100_C40,15000,32768,130,0,0,0,0);
data100A(07,:)=phase2(CF3_15_90_100_C50,15000,32768,130,0,0,0,0);
data100A(08,:)=phase2(CF3_15_90_100_C55,15000,32768,130,0,0,0,0);
data100A(09,:)=phase2(CF3_15_90_100_C60,15000,32768,130,0,0,0,0);
data100A(10,:)=phase2(CF3_15_90_100_C70,15000,32768,130,0,0,0,0);
data100A(11,:)=phase2(CF3_15_90_100_C80,15000,32768,130,0,0,0,0);
data100A(12,:)=phase2(CF3_15_90_100_C90,15000,32768,130,0,0,0,0);

data420A(01,:)=phase2(CF3_15_90_420_C00,15000,32768,130,0,0,0,0);
data420A(02,:)=phase2(CF3_15_90_420_C10,15000,32768,130,0,0,0,0);
data420A(03,:)=phase2(CF3_15_90_420_C20,15000,32768,130,0,0,0,0);
data420A(04,:)=phase2(CF3_15_90_420_C25,15000,32768,130,0,0,0,0);
data420A(05,:)=phase2(CF3_15_90_420_C30,15000,32768,130,0,0,0,0);
data420A(06,:)=phase2(CF3_15_90_420_C40,15000,32768,130,0,0,0,0);
data420A(07,:)=phase2(CF3_15_90_420_C50,15000,32768,130,0,0,0,0);
data420A(08,:)=phase2(CF3_15_90_420_C55,15000,32768,130,0,0,0,0);
data420A(09,:)=phase2(CF3_15_90_420_C60,15000,32768,130,0,0,0,0);
data420A(10,:)=phase2(CF3_15_90_420_C70,15000,32768,130,0,0,0,0);
data420A(11,:)=phase2(CF3_15_90_420_C80,15000,32768,130,0,0,0,0);
data420A(12,:)=phase2(CF3_15_90_420_C90,15000,32768,130,0,0,0,0);

data312A(01,:)=phase2(CF3_15_90_312_C00,15000,32768,130,0,0,0,0);
data312A(02,:)=phase2(CF3_15_90_312_C10,15000,32768,130,0,0,0,0);
data312A(03,:)=phase2(CF3_15_90_312_C20,15000,32768,130,0,0,0,0);
data312A(04,:)=phase2(CF3_15_90_312_C25,15000,32768,130,0,0,0,0);
data312A(05,:)=phase2(CF3_15_90_312_C30,15000,32768,130,0,0,0,0);
data312A(06,:)=phase2(CF3_15_90_312_C40,15000,32768,130,0,0,0,0);
data312A(07,:)=phase2(CF3_15_90_312_C50,15000,32768,130,0,0,0,0);
data312A(08,:)=phase2(CF3_15_90_312_C55,15000,32768,130,0,0,0,0);
data312A(09,:)=phase2(CF3_15_90_312_C60,15000,32768,130,0,0,0,0);
data312A(10,:)=phase2(CF3_15_90_312_C70,15000,32768,130,0,0,0,0);
data312A(11,:)=phase2(CF3_15_90_312_C80,15000,32768,130,0,0,0,0);
data312A(12,:)=phase2(CF3_15_90_312_C90,15000,32768,130,0,0,0,0);

data160A(01,:)=phase2(CF3_15_90_160_C00,15000,32768,130,0,0,0,0);
data160A(02,:)=phase2(CF3_15_90_160_C10,15000,32768,130,0,0,0,0);
data160A(03,:)=phase2(CF3_15_90_160_C20,15000,32768,130,0,0,0,0);
data160A(04,:)=phase2(CF3_15_90_160_C25,15000,32768,130,0,0,0,0);
data160A(05,:)=phase2(CF3_15_90_160_C30,15000,32768,130,0,0,0,0);
data160A(06,:)=phase2(CF3_15_90_160_C40,15000,32768,130,0,0,0,0);
data160A(07,:)=phase2(CF3_15_90_160_C50,15000,32768,130,0,0,0,0);
data160A(08,:)=phase2(CF3_15_90_160_C55,15000,32768,130,0,0,0,0);
data160A(09,:)=phase2(CF3_15_90_160_C60,15000,32768,130,0,0,0,0);
data160A(10,:)=phase2(CF3_15_90_160_C70,15000,32768,130,0,0,0,0);
data160A(11,:)=phase2(CF3_15_90_160_C80,15000,32768,130,0,0,0,0);
data160A(12,:)=phase2(CF3_15_90_160_C90,15000,32768,130,0,0,0,0);

data670A(01,:)=phase2(CF3_15_90_670_C00,15000,32768,130,0,0,0,0);
data670A(02,:)=phase2(CF3_15_90_670_C10,15000,32768,130,0,0,0,0);
data670A(03,:)=phase2(CF3_15_90_670_C20,15000,32768,130,0,0,0,0);
data670A(04,:)=phase2(CF3_15_90_670_C25,15000,32768,130,0,0,0,0);
data670A(05,:)=phase2(CF3_15_90_670_C30,15000,32768,130,0,0,0,0);
data670A(06,:)=phase2(CF3_15_90_670_C40,15000,32768,130,0,0,0,0);
data670A(07,:)=phase2(CF3_15_90_670_C50,15000,32768,130,0,0,0,0);
data670A(08,:)=phase2(CF3_15_90_670_C55,15000,32768,130,0,0,0,0);
data670A(09,:)=phase2(CF3_15_90_670_C60,15000,32768,130,0,0,0,0);
data670A(10,:)=phase2(CF3_15_90_670_C70,15000,32768,130,0,0,0,0);
data670A(11,:)=phase2(CF3_15_90_670_C80,15000,32768,130,0,0,0,0);
data670A(12,:)=phase2(CF3_15_90_670_C90,15000,32768,130,0,0,0,0);

theta=[0 10 20 25 30 40 50 55 60 70 80 90]; 

xax=((0:32767)/32768-0.5)*15;
fig = figure(1);
for k=1:12
  subplot(3,4,k)
  plot(xax,real(data160A(k,:)));
  axis([-2 2 0 25000])
  set(gca,'FontSize',14)
  xlabel('\nu [kHz]')
  line = sprintf('\\theta = %02d\\circ',theta(k));
  title(line)
  set(gca, 'XDir','reverse')

end
orient('landscape')
print(fig, '-fillpage', '-dpdf', 'new_figure_C_160.pdf')

xax=((0:32767)/32768-0.5)*15;
fig = figure(2);
for k=1:12
  subplot(3,4,k)
  plot(xax,real(data670A(k,:)));
  axis([-2 2 0 25000])
  set(gca,'FontSize',14)
  xlabel('\nu [kHz]')
  line = sprintf('\\theta = %02d\\circ',theta(k));
  title(line)
  set(gca, 'XDir','reverse')

end
orient('landscape')
print(fig, '-fillpage', '-dpdf', 'new_figure_C_670.pdf')

fig = figure(3);
clf
plot(xax,data160A(12,:)+5e3*0,'r',xax,data670A(12,:)+5e3*0,'k')
hold on
plot(xax,data160A(08,:)+5e3*1,'r',xax,data670A(08,:)+5e3*1,'k')
plot(xax,data160A(06,:)+5e3*2,'r',xax,data670A(06,:)+5e3*2,'k')
plot(xax,data160A(05,:)+5e3*3,'r',xax,data670A(05,:)+5e3*3,'k')
plot(xax,data160A(04,:)+5e3*4,'r',xax,data670A(04,:)+5e3*4,'k')
plot(xax,data160A(03,:)+5e3*5,'r',xax,data670A(03,:)+5e3*5,'k')
hold off
axis([-2 2 0 3.5e4])
set(gca, 'XDir','reverse')

orient('tall')
print(fig, '-dpdf', '-fillpage', 'new_figure_paperlike2.pdf')


fig = figure(4);
clf
plot(xax,data160A(12,:)+5e3*0,'r',xax,data420A(12,:)+5e3*0,'k')
hold on
plot(xax,data160A(08,:)+5e3*1,'r',xax,data420A(08,:)+5e3*1,'k')
plot(xax,data160A(06,:)+5e3*2,'r',xax,data420A(06,:)+5e3*2,'k')
plot(xax,data160A(05,:)+5e3*3,'r',xax,data420A(05,:)+5e3*3,'k')
plot(xax,data160A(04,:)+5e3*4,'r',xax,data420A(04,:)+5e3*4,'k')
plot(xax,data160A(03,:)+5e3*5,'r',xax,data420A(03,:)+5e3*5,'k')
hold off
axis([-2 2 0 3.5e4])
set(gca, 'XDir','reverse')


orient('tall')
print(fig, '-dpdf', '-fillpage', 'new_figure_paperlike.pdf')

% Plot the three k_ex values representing rac and S (S is 3 or 4 times big)
fig = figure(5);
clf
%theta=[0 10 20 25 30 40 50 55 60 70 80 90]; 
%90 degrees
plot(xax,data160A(12,:)+5e3*0,'r',xax,data420A(12,:)+5e3*0,'k', xax, data670A(12,:)+5e3*0, 'b')
hold on
% 55 degrees
plot(xax,data160A(08,:)+5e3*1,'r',xax,data420A(08,:)+5e3*1,'k', xax, data670A(08,:)+5e3*1, 'b')
% 40 degrees
plot(xax,data160A(06,:)+5e3*2,'r',xax,data420A(06,:)+5e3*2,'k', xax, data670A(06,:)+5e3*2, 'b')
% 30 degrees
plot(xax,data160A(05,:)+5e3*3,'r',xax,data420A(05,:)+5e3*3,'k', xax, data670A(05,:)+5e3*3, 'b')
% 25 degrees
plot(xax,data160A(04,:)+5e3*4,'r',xax,data420A(04,:)+5e3*4,'k', xax, data670A(04,:)+5e3*4, 'b')
% 20 degrees
plot(xax,data160A(03,:)+5e3*5,'r',xax,data420A(03,:)+5e3*5,'k', xax, data670A(03,:)+5e3*5, 'b')
hold off
axis([-2 2 0 3.5e4])
set(gca, 'XDir','reverse')


orient('tall')
print(fig, '-dpdf', '-fillpage', 'plot_three_kex_vals_160_420_670.pdf')


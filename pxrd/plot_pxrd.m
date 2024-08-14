% plot pxrd patterns and take difference

simulated_s_tfla = readtable("666327_S_TFLA.txt");
experimental_s_tfla = readtable("D8-PDF-000471_exported.txt");

figure(1); hold on;
xsim = simulated_s_tfla{:,1} -0.2;
ysim = simulated_s_tfla{:,2};
plot(xsim, ysim)

xexp = experimental_s_tfla{:,1};
yexp = experimental_s_tfla{:,2};
plot(xexp, yexp)

xlim([5 25])
hold off;

% shortysim = ysim(1:4052);
% diff = yexp - shortysim;
% figure(1111)
% plot(xexp, diff)
% 
% 
% 
% xsim_Q = ((4*pi)/0.55941)*(sin(pi.*xsim)/360);
% xexp_Q = ((4*pi)/0.55941)*(sin(pi.*xexp)/360);
% 
% figure(2); hold on;
% plot(xsim_Q, ysim)
% plot(xexp_Q, yexp)
% 
% hold off;

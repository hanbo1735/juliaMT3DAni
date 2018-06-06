%%
close all; fclose all; clear; clc

%%
RespSpec = {'*.resp','Resp Files(*.resp)';'*.*','All Files(*.*)'};

[filen, path] = uigetfile(RespSpec,'Open the analytic response file:');
respfile = fullfile(path, filen);
resp0 = readMT1DAniResp(respfile);


[filen, path] = uigetfile(RespSpec, 'Open the numerical response file:');
respfile = fullfile(path, filen);
[datInfo, resp1] = readMT3DResp(respfile);

%% select data
periods = 1 ./ resp0.freqs;

rhoxy0 =  resp0.appRho(:, 3:3);
phsxy0 = -resp0.appRho(:, 4:4);
rhoyx0 =  resp0.appRho(:, 5:5);
phsyx0 = -resp0.appRho(:, 6:6);

siteNo = round(size(datInfo.rxLoc, 1)/2);
subRxID = find(datInfo.rxID==siteNo);

rhoxy1 =  resp1(subRxID, 3:3);
phsxy1 =  resp1(subRxID, 4:4);
rhoyx1 =  resp1(subRxID, 5:5);
phsyx1 =  resp1(subRxID, 6:6);

%% plot
% apparent resistivity
figure(1);

p1 = loglog(periods, rhoxy0,'-b','linewidth',1);
hold on;
p2 = loglog(periods, rhoxy1,'ob','MarkerFaceColor','b','MarkerSize',5,'linewidth',2,'linestyle','none');

p3 = loglog(periods, rhoyx0,'-r','linewidth',1);
p4 = loglog(periods, rhoyx1,'or','MarkerFaceColor','r','MarkerSize',5,'linewidth',2,'linestyle','none');

labelx = 'Period (s)';
labely_amp = 'App.Res (Ohm)';

xlabel(labelx,'fontsize',15);
ylabel(labely_amp,'fontsize',15);

h_leg = legend([p1 p3],'XY','YX');

set(gca,'fontsize',15, 'linewidth',1.5);
set(h_leg,'fontsize',14, 'linewidth',1);
% legend('boxoff');

set(gca,'ticklength',[0.02 0.02]);
set(gca,'xminortick','on');
% grid on;


%%
% apparent resistivity error
figure(3);
rhoxy_e = abs(rhoxy1-rhoxy0) ./ abs(rhoxy0) * 100;
rhoyx_e = abs(rhoyx1-rhoyx0) ./ abs(rhoyx0) * 100;

% p1 = semilogx(periods, rhoxy_e,'sb','MarkerFaceColor','b','MarkerSize',5,'linewidth',2,'linestyle','none');
p1 = semilogx(periods, rhoxy_e,'-b','linewidth',2);
hold on;
% p2 = semilogx(periods, rhoyx_e,'sr','MarkerFaceColor','r','MarkerSize',5,'linewidth',2,'linestyle','none');
p2 = semilogx(periods, rhoyx_e,'-r','linewidth',2);


% p2 = loglog(freqs, rhoxy1,'ob','MarkerFaceColor','b','MarkerSize',5,'linewidth',2,'linestyle','none');
% 
% p3 = loglog(freqs, rhoyx0,'-r','linewidth',1);
% p4 = loglog(freqs, rhoyx1,'or','MarkerFaceColor','r','MarkerSize',5,'linewidth',2,'linestyle','none');

labelx = 'Period (s)';
labely_amp = 'App.Res relative error (%)';

xlabel(labelx,'fontsize',15);
ylabel(labely_amp,'fontsize',15);

% h_leg = legend([p1 p2 p3],'Ey','Bx','Ez');

set(gca,'fontsize',15, 'linewidth',1.5);
% set(h_leg,'fontsize',13, 'linewidth',1);
legend('boxoff');
% set(gca,'xminorgrid','on');


%%
% impedance phase
figure(2);
p1 = semilogx(periods, phsxy0,'-b','linewidth',1);
hold on;
p2 = semilogx(periods, phsxy1,'ob','MarkerFaceColor','b','MarkerSize',5,'linewidth',2,'linestyle','none');

p3 = semilogx(periods, phsyx0,'-r','linewidth',1);
p4 = semilogx(periods, phsyx1,'or','MarkerFaceColor','r','MarkerSize',5,'linewidth',2,'linestyle','none');

labelx = 'Period (s)';
labely_pha = 'Phase (\circ)';

xlabel(labelx,'fontsize',15);
ylabel(labely_pha,'fontsize',15);

set(gca,'fontsize',15, 'linewidth',1.5);

% h_leg = legend([p1 p3],'\it\phi_x_y','\it\phi_y_x');
% set(h_leg,'fontsize',14, 'linewidth',1);
% legend('boxoff');


%%
% phase error
figure(4);
phsxy_e = phsxy1 - phsxy0;
phsyx_e = phsyx1 - phsyx0;
% p1 = semilogx(periods, phsxy_e,'sb','MarkerFaceColor','b','MarkerSize',5,'linewidth',2,'linestyle','none');
p1 = semilogx(periods, phsxy_e,'-b','linewidth',2);
hold on;
% p2 = semilogx(periods, phsyx_e,'sr','MarkerFaceColor','r','MarkerSize',5,'linewidth',2,'linestyle','none');
p2 = semilogx(periods, phsyx_e,'-r','linewidth',2);

labelx = 'Period (s)';
labely_pha = 'Phase difference (\circ)';

xlabel(labelx,'fontsize',15);
ylabel(labely_pha,'fontsize',15);

% h_leg = legend([p1 p2 p3],'Ey','Bx','Ez');

set(gca,'fontsize',15, 'linewidth',1.5);
% set(h_leg,'fontsize',13, 'linewidth',1);
legend('boxoff');
% set(gca,'xminorgrid','on');




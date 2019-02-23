% function relatepf2pfix1
close all;
colors=[0, 0, 128;
    0,191,255;
    0, 128, 128;
    0, 128, 0;
%     255, 255, 0;
    128, 128, 0;
    255, 165, 0;
    255,69,0;
    205, 92,92;
    255, 0, 0;
    128, 0, 0]./255;
Markers=['*','.','o','>','p','x','<'];
Pf1s=[0.0005:0.0001:0.0009]';
Pfix1=0.01:0.01:0.99;
Pf2s=zeros(length(Pf1s),length(Pfix1));
fontsize=15;
figure(1);

Pftt=0.001;
for i=1:length(Pf1s)
    Pf2s(i,:)=(Pftt-Pf1s(i))./(1-Pfix1);
%     plot(Pfix1,100*Pf2s(i,:),'-','Color',colors(i,:),'Marker',Markers(i)); hold on;
        semilogy(Pfix1,Pf2s(i,:),'-','Color',colors(i,:),'Marker',Markers(i)); hold on;

end
legend('P_f_1=0.0005','P_f_1=0.0006','P_f_1=0.0007','P_f_1=0.0008','P_f_1=0.0009');
xlabel('P_f_i_x_1','FontSize',fontsize);
ylabel('P_f_2','FontSize',fontsize);
set(gca,'FontSize',fontsize);
% ,'Color',colors(i,:),'Marker',Markers(i)
figure(2);
for i=1:length(Pfix1)
    Pf2s(:,i)=(Pftt-Pf1s)/(1-Pfix1(i));
    semilogy(Pf1s,Pf2s(:,i),'-*'); hold on;

end

N = [2 3 4 5 6 7];
distributed = [420 3.718000e+02 3.382000e+02 3.091000e+02 2.838000e+02 2.595000e+02];
classicalpamArray1 = [406 3.573000e+02 3.173000e+02 2.875000e+02 2.617000e+02 238];
fastpamArray1 = [406 3.581000e+02 3.156000e+02 2.861000e+02 2.617000e+02 2.381000e+02];
newpamArray1 = [406 355 317 287 263 238];
%distributed = [9.037400e+02 7.859800e+02 7.101200e+02 6.312600e+02 5.776500e+02 5.341100e+02];
%classicalpamArray1 = [816 698 608 540 4.933400e+02 4.577100e+02];
%fastpamArray1 = [8.284200e+02 7.124900e+02 6.232700e+02 5.581200e+02 5.082600e+02 4.712900e+02];
%newpamArray1 = [816 698 608 540 493 458];
XTickString1 = cell(1,length(N));
for i = 1 : length(N)
    temp1 = num2str(N(1,i));
    XTickString1(1,i) = {temp1};

end
figure('units','inches')
pos = get(gcf,'pos');
set(gcf,'pos',[pos(1) pos(2) 9 4])
XTickString11 = categorical(XTickString1);
XTickString1 = reordercats(XTickString11,XTickString1);
b = bar(XTickString1, [distributed;classicalpamArray1;fastpamArray1; newpamArray1], 0.8);
leg = legend('DHMCA', 'PAMwithkmeans++','fast\_kmedoids', 'improved\_kmedoids');
leg.Location = 'northeast';
set(leg,'FontSize',10)
ylim([0 1200])
ax = gca;
ax.FontSize = 11;
xlabel('number of clusters', 'FontSize',14)
ylabel('performance', 'FontSize',14)
%print -depsc2 comparisonDitributedWaterNetwork


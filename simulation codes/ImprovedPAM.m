function [idx,C,TD,ite] = ImprovedPAM(X,k)
% 输入：X为聚类数据矩阵（每行为一个点），k为聚类中心个数
% 输出：idx为每个点的聚类标签，C为聚类中心，TD为Total deviation loss，ite为迭代次数
% 参考文献：Schubert, E., & Rousseeuw, P. J. (2021). Fast and eager k-medoids clustering: O (k) runtime improvement of the PAM, CLARA, and CLARANS algorithms. Information Systems, 101, 101804.

%%% 数据点差异性计算
%D = squareform(pdist(X,'euclidean')); index = 1:length(D);
D = X;index = 1:length(D);
C = nan*ones(k,1);
%%% PAM BUILD: Find initial cluster centers
[TD,C(1)] = min(sum(D,2)-diag(D));
dno = D(:,C(1));
for i = 2:k
    delta_TD_star = inf; x_star = nan;
    XC = setdiff(index,C(1:i-1));
    for xc = 1:length(XC)
        delta_TD = 0;
        for xo = 1:length(XC)
            delta = D(XC(xo),XC(xc))-dno(XC(xo));
            if delta < 0
                delta_TD = delta_TD + delta;
            end
        end
        if delta_TD < delta_TD_star
            delta_TD_star = delta_TD; x_star = XC(xc);
        end
    end
    TD = TD + delta_TD_star; C(i) = x_star;
    XO = setdiff(XC,C(i));
    dno(XO) = min([dno(XO)';D(XO,C(i))']);
end
%%% FastPAM1: Improved SWAP algorithm
XC = setdiff(index,C); [dnearXO, nearXO] = min(D(:,C),[],2); 
TempD = D(:,C);
for i = 1:length(index)
    TempD(i,nearXO(i)) = inf;
end
dsecondXO = min(TempD,[],2); 
delta_TD_C = nan*ones(k,1);
ite = 0;
% disp(['第' num2str(ite) '次迭代TD：' num2str(TD)]);
% disp('中心：');
% disp(C);
while true
    ite = ite + 1;
    for i = 1:k
        nearXOC = find(nearXO==i);
        delta_TD_C(i) = sum(dsecondXO(nearXOC)-dnearXO(nearXOC));
    end
    delta_TD_star = 0; m_star = nan; x_star = nan;
    for xc = 1:length(XC)
        delta_TD = delta_TD_C;
        delta_TD_xc = 0;
        for xo = 1:length(index)
            doj = D(xo,XC(xc));
            if doj < dnearXO(xo)
                delta_TD_xc = delta_TD_xc + doj - dnearXO(xo);
                CnearXO = find(C == C(nearXO(xo)));
                delta_TD(CnearXO) = delta_TD(CnearXO) + dnearXO(xo) - dsecondXO(xo);
            elseif doj < dsecondXO(xo)
                CnearXO = find(C == C(nearXO(xo)));
                delta_TD(CnearXO) = delta_TD(CnearXO) + doj - dsecondXO(xo);
            else
            end
        end
        [~,i] = min(delta_TD);
        delta_TD(i) = delta_TD(i) + delta_TD_xc;
        if delta_TD(i) < delta_TD_star
            delta_TD_star = delta_TD(i);
            m_star = C(i);
            x_star = XC(xc);
        end
    end
    if delta_TD_star >= 0
        break;
    end
    C(C==m_star) = x_star;
    XC = setdiff(index,C); [dnearXO, nearXO] = min(D(:,C),[],2); 
    TempD = D(:,C);
    for i = 1:length(index)
        TempD(i,nearXO(i)) = inf;
    end
    dsecondXO = min(TempD,[],2); 
    TD = TD + delta_TD_star;
end
idx = C(nearXO);
end


function [out] = spc_summary(result, MyEvn)
% summary function for SPC.m

useiters =(MyEvn.burnin + 1):MyEvn.tot; %second half for posterior inference
nSample = numel(useiters); 
ntot = nSample*MyEvn.nchain; %nIter - burnin;
N = size(MyEvn.W, 1);  p = size(result{1}.BETA{1},1);
Nus = nan(ntot, N);  AllBeta = zeros(N, p);
WeightMatr = zeros(N);  Allds = nan(1, ntot);  Allcenter = cell(1,ntot);
Allcluster = zeros(ntot, N);  Allbetar = cell(1,ntot);
for ch = 1:MyEvn.nchain
    cind = (ch-1)*nSample + (1:nSample); 
    Allds(cind) = result{ch}.Num(useiters);
    Allcluster(cind,:) = result{ch}.Cluster(useiters, :);
    Nus(cind, :) = result{ch}.Nu(useiters,:);
    for i = 1:nSample
        labs = zeros(1,N);
        Allcenter{(ch-1)*nSample + i} = result{ch}.Center{useiters(i)};
        Allbetar{(ch-1)*nSample + i} = result{ch}.BETA{useiters(i)};
        for j = 1:length(result{ch}.Center{useiters(i)})
            inds = find(result{ch}.Cluster(useiters(i),:) == result{ch}.Center{useiters(i)}(j));
            labs(inds) = j;
            WeightMatr(inds, inds) = WeightMatr(inds, inds) + 1;
        end
        for j = 1:p
            AllBeta(:,j) = AllBeta(:,j) + result{ch}.BETA{useiters(i)}(j,labs)';
        end
    end
end
WeightMatr = WeightMatr/ntot;
for i = 1:N;  WeightMatr(i,i) = 1;  end
count = histc(Allds, 1:max(Allds));   Kmode = find(count == max(count));
ClusterOut = JWClus(WeightMatr, Kmode);

mat_REGC = nan(N,3);
for i = 1:N
    tmp = prctile(Nus(:,i), [5 95]);
    mat_REGC(i, :) = [mean(Nus(:,i)), tmp];
end

%-------------------------- get beta_r
BETAmatr = [];
index=find(Allds==Kmode); nindex=length(index);
for i = 1:nindex
    shind = index(i);
    BETA = Allbetar{shind};
    Center = Allcenter{shind};
    Cluster = Allcluster(shind,:);
    [map, common] = ClusMap(Cluster,Center,ClusterOut);
    for j = 1:length(Center)
        conindex = find(map==j);
        BETAmatr(i,j,:)=BETA(:,conindex);
    end
end
% K by p summary
BetaR_mean = squeeze(mean(BETAmatr,1));
BetaR_lb = squeeze(quantile(BETAmatr,0.05,1));
BetaR_ub = squeeze(quantile(BETAmatr,0.95,1));

AllBeta = AllBeta./ntot;
Nus_obs = log((MyEvn.Observed'+0.5*(MyEvn.Observed'==0))./MyEvn.Expected'); %not accurate
err1 = repmat(Nus_obs, [ntot,1]) - Nus;  %error terms REGC

% moran's I for posterior samples of residuals
Spatind = nan(1, ntot);    for i = 1:ntot;  Spatind(i) = SpatStat(err1(i,:)',MyEvn.W,1,1);   end

% save output
out.labs = ClusterOut; %clustering label
out.K = Kmode; %posterior mode of number of clusters
% K by p summary of cluster-wise regression coefficients
out.BetaR_mean = BetaR_mean; 
out.BetaR_lb = BetaR_lb; 
out.BetaR_ub = BetaR_ub; 
out.Spatind = Spatind; %moran's I for posterior samples of residuals
end

function [val, Z] = SpatStat(X, W, pop, ind)
% [val] = SpatStat(X,W,pop,ind)
% function to calculate the spatial dependence statistics
% X: the n by 1 vector of observed data nad pop its the population  
% W: n by n neighborhood matrix
% ind = 1:moran's I   2:Geary's C   3:Oden's Ipop   4:Tango
n=size(W,1); pop = sum(pop);
meanX=mean(X);
switch ind
    case{1} %moran's I
        val = sum(sum(W.*((X-meanX)*(X-meanX)')))/sum(diag((X-meanX)*(X-meanX)'));
        val = val*(n/sum(sum(W)));
        Ms = -1/(n - 1);
        S1 = .5*sum(sum((W+W').^2));
        S2 = sum((sum(W,1) + sum(W,2)').^2);
        S3 = sum((X-meanX).^4)/n/((sum((X-meanX).^2)/n)^2);
        S4 = (n^2-3*n+3)*S1 - n*S2 + 3*(sum(sum(W)))^2;
        S5 = S1 - 2*n*S1 + 6*(sum(sum(W)))^2;
        Vars = (n*S4 - S3*S5)/((n-1)*(n-2)*(n-3)*(sum(sum(W)))^2);
        Z = (val - Ms)/sqrt(Vars);
               
    case{2} %Geary's C (3.13)
        % values: 0 to 2 with 1: no correlation, 2:strong negative relationship, 0: strong positive relationship
        val = sum(sum(W.*(repmat(X,1,n)-repmat(X',n,1)).^2))/sum(diag((X-meanX)*(X-meanX)'));
        val = val*((n-1)/2/(sum(sum(W)) - sum(diag(W))));
    case{3} % Oden's Ipop statistics (3.15)
        p = 1/n*ones(n,1); % assume expected value is uniformly distributed
        r = (X./pop)/sum(X./pop);
        % pop is the total population
        bbar = sum(X)/pop; 
        A = sum(sum(W.*(p*p'))); B = sum(diag(W).*p);
        S0 = pop^2*A-pop*B;
        val = (sum(X))^2*sum(sum(W.*((r-p)*(r-p)'))) - sum(X)*(1-2*bbar)*sum(diag(W).*r) - ...
            sum(X)*bbar*sum(diag(W).*p);
        val = val/(S0*bbar*(1-bbar));
    case{4} % Tango (3.22)
        p = 1/n*ones(n,1); % assume expeced is uniformly distributed
        r = (X./pop)/sum(X./pop);
        val = sum(sum(W.*((r-p)*(r-p)')));
    %case{5} % Getis and Ord's global statistics - G(d) (3.32), need specify distance d
    %    % the sum is over pairs of regions that have distance measure <
    %    % threshold
    %    val = sum(sum(W.*(X*X')))/sum(sum(X*X'));
        
    otherwise
        display('index not found.')
end
end

function [IDX, C] = JWClus(affinity, K)
for i=1:size(affinity,1);  D(i,i) = sum(affinity(i,:));  end
for i=1:size(affinity,1)
    for j=1:size(affinity,2)
        NL1(i,j) = affinity(i,j) / (sqrt(D(i,i)) * sqrt(D(j,j)));  
    end
end
[eigVectors, eigValues] = eig(NL1);
k = K;
nEigVec = eigVectors(:,(size(eigVectors,1)-(k-1)): size(eigVectors,1));
for i=1:size(nEigVec,1)
    n = sqrt(sum(nEigVec(i,:).^2));    
    U(i,:) = nEigVec(i,:) ./ n; 
end
[IDX, C] = kmeans(U,K);
end

function [map, common] = ClusMap( ResultOut,Center,ClusterOut )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
K=length(Center);
ComMat=NaN(K,K);
for i=1:K
    for j=1:K
        ResIndex=find(ResultOut==Center(i));
        ClsIndex=find(ClusterOut==j);
        ComMat(i,j)=length(intersect(ResIndex,ClsIndex));
    end
end
map=NaN(1,K);
common=NaN(1,K);
p=max(ComMat(:));
[x, y]=find(ComMat==p);
t = length(x);
if t==1
    map(x)=y;
    common(x)=p;
else 
    sam=randsample(t,1);
    map(x(sam))=y(sam);
    common(x(sam))=p;
end
for l=2:(K-1)
    rowexist=find(map>0);
    rowIndex=setdiff(1:K,rowexist);
    colexist=map(~isnan(map));
    colIndex=setdiff(1:K,colexist);
    subMat=ComMat(rowIndex,colIndex);
    q=max(subMat(:));
    [x y]=find(subMat==q);
    t=length(x);
if t==1
    map(rowIndex(x))=colIndex(y);
    common(rowIndex(x))=q;
else 
    sam=randsample(t,1);
    a=x(sam);
    b=y(sam);
    map(rowIndex(a))=colIndex(b);
    common(rowIndex(a))=q;
end
end
common(find(isnan(map)==1))=ComMat(find(isnan(map)==1),setdiff(1:K,map(~isnan(map))));
map(isnan(map)==1)=setdiff(1:K,map(~isnan(map)));
end






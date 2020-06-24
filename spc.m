function [out] = spc(ev)
% main function for spatial clustering

ev.DistMatr = DistMatrix(ev.Latitude, ev.Longitude); %compute distance
iniK = 5;  %initial number of clusters
d = 0.01;  %penalty parameter in the prior cluster model
MyIni = initPara(d, iniK, ev);
MySim = MySimFun(ev);

bloop = 1;
New = MyIni;
tic
starttime = cputime;
while bloop <= ev.tot;
     if ev.verbose;  fprintf('%6d', bloop);  if(~mod(bloop, 10));  fprintf('\n');  end;  end
    Old = New;
    u = rand(1);
    if Old.Num ~= 1
        if u <= 0.4
            [New,B,C,D] = BirthMove(Old,ev, MySim, bloop);
            bloop = C;
            out.alpha(bloop) = B;
            out.multi(bloop) = D;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 1;
            num = New.Num;
            indic = 1;
            Old = New;
            [New,B,C,D] = HeightMove(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = B;
            out.multi(bloop) = D;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 3;
            num = New.Num;
            indic = 3;
            Old = New;
            [New,C] = HyperMove(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = 1;
            out.multi(bloop) = 1;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 4;
            num = New.Num;
            indic = 4;
            Old = New;
            [New,C] = updateD(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = 1;
            out.multi(bloop) = 1;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 6;
            num = New.Num;
            indic = 6;
        elseif  u>0.4 && u<=0.8
            [New,B,C,D] = DeathMove(Old,ev,MySim,bloop);
            bloop=C;
            out.alpha(bloop) = B;
            out.multi(bloop) = D;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 2;
            num = New.Num;
            indic = 2;
            Old = New;
            [New,B,C,D] = HeightMove(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = B;
            out.multi(bloop) = D;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 3;
            num = New.Num;
            indic = 3;
            Old=New;
            [New,C] = HyperMove(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = 1;
            out.multi(bloop) = 1;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 4;
            num = New.Num;
            indic = 4;
            Old = New;
            [New,C] = updateD(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = 1;
            out.multi(bloop) = 1;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 6;
            num = New.Num;
            indic = 6;
        else [New,B,C,D] = ShiftMove(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = B;
            out.multi(bloop) = D;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 5;
            num = New.Num;
            indic = 5;
        end
    else
        if u <= 0.8
            [New,B,C,D] = BirthMove(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = B;
            out.multi(bloop) = D;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 1;
            num = New.Num;
            indic = 1;
            Old = New;
            [New,B,C,D] = HeightMove(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = B;
            out.multi(bloop) = D;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 3;
            num = New.Num;
            indic = 3;
            Old = New;
            [New,C] = HyperMove(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = 1;
            out.multi(bloop) = 1;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 4;
            num = New.Num;
            indic = 4;
            Old = New;
            [New,C] = updateD(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = 1;
            out.multi(bloop) = 1;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 6;
            num = New.Num;
            indic = 6;
        else [New,B,C,D] = ShiftMove(Old,ev,MySim,bloop);
            bloop = C;
            out.alpha(bloop) = B;
            out.multi(bloop) = D;
            out.Nu(bloop,:) = New.Nu;
            out.Cluster(bloop,:) = New.Cluster;
            out.Num(bloop) = New.Num;
            out.sigma2(bloop) = New.sigma2;
            out.BETA{bloop} = New.BETA;
            out.Sigma2{bloop} = New.Sigma2;
            out.Center{bloop} = New.Center;
            out.d(bloop) = New.d;
            out.flag(bloop) = 5;
            num = New.Num;
            indic = 5;
        end
    end
end
endtime = cputime;
out.usedtime = toc;
out.duration = endtime - starttime;
CPUtime = out.usedtime/3600;
fprintf('\n%d iterations are done with elapsed time %.2f hours.\n', ev.tot, CPUtime)
% save(strcat('realresult',num2str(ncas),'_',num2str(ch),'.mat'),'result', 'ev', 'MyIni')
end

function [ ini ] = initPara(d, iniK, MyEvn)
Coordinate = [MyEvn.Latitude MyEvn.Longitude];
IDX = kmeans(Coordinate, iniK);
minsize = size(MyEvn.X,2);
while minsize <= min(size(MyEvn.X, 2), MyEvn.nrmin)
    iniCenter = NaN(1,iniK);
    for i = 1:iniK;  iniCenter(i)=randsample(find(IDX==i),1);  end
    ini.Center = MyEvn.DistIndex(iniCenter);
    ini.Cluster = ClusterGen(MyEvn.DistMatr, ini.Center, MyEvn.DistIndex);
    Clustersize = NaN(1,iniK);
    for i = 1:iniK;   Clustersize(i) = length(find(ini.Cluster==ini.Center(i)));   end
    minsize = min(Clustersize);
end
ini.Num = iniK;
ini.Nu = log((MyEvn.Observed'+0.5*(MyEvn.Observed'==0))./MyEvn.Expected');
ini.BETA = NaN(size(MyEvn.X,2),iniK);
ini.Sigma2 = NaN(1,iniK);
sigma2 = NaN(1,iniK);
for i = 1:iniK
    ind = find(ini.Cluster==ini.Center(i));
    [b,bint,r,rint,stats] = regress(ini.Nu(ind)',MyEvn.X(ind,:));
    ini.BETA(:,i) = b;
    sigma2(i) = stats(4);
end
ini.sigma2 = 1/2*min(sigma2);
for i = 1:iniK;  ini.Sigma2(i)=sigma2(i)-ini.sigma2;  end
ini.d = d;
end

function [ MySim ] = MySimFun(MyEvn)
MySim.Obs = MyEvn.Observed';
MySim.obsNu = log((MySim.Obs+0.5)'./MyEvn.Expected);
end

function [ New,alpha,bloop,multi,ARatio,PRatio ] = BirthMove( Old,MyEvn,MySim,bloop )
Update.Num = Old.Num+1;
newCenter = randsample(setdiff(MyEvn.DistIndex,Old.Center),1);
Update.Center = NaN(1,Update.Num);
ind = randsample(Update.Num,1);
Update.Center(ind) = newCenter;
Update.Center(setdiff(1:Update.Num,ind)) = Old.Center;
Update.Cluster = ClusterGen(MyEvn.DistMatr,Update.Center,MyEvn.DistIndex);
count = histc(Update.Cluster, sort(unique(Update.Cluster)));
New = Old; alpha = -99; multi = -99; ARatio = 0; PRatio = 0;
if min(count) >= MyEvn.nrmin
    newIndex = find(Update.Cluster==newCenter);
    %Update from original proposal
    [UpdateSigma2, UpdateBETA]=UpdateNew(Update.Cluster,newCenter,Old,MyEvn,MySim);
    %Update from prior proposal
    %[UpdateSigma2 UpdateBETA]=UpdateNew2(Update.Cluster,newCenter,Old,MyEvn,MySim);
    Update.Nu=Old.Nu;
    Update.Sigma2=NaN(1,Update.Num);
    Update.Sigma2(ind)=UpdateSigma2;
    Update.Sigma2(setdiff(1:Update.Num,ind))=Old.Sigma2;
    Update.BETA=NaN(size(Old.BETA,1),Update.Num);
    Update.BETA(:,ind)=UpdateBETA;
    Update.BETA(:,setdiff(1:Update.Num,ind))=Old.BETA;
    Update.sigma2=Old.sigma2;
    Update.d=Old.d;
    Like=1;
    ratio=NuDensity(Update.Nu,Update.Center,Update.Cluster,Update.BETA,Update.Sigma2,Update.sigma2,MyEvn.X)/NuDensity(Old.Nu,Old.Center,Old.Cluster,Old.BETA,Old.Sigma2,Old.sigma2,MyEvn.X);
    ARatio=(1-Update.d)/(MyEvn.n-Old.Num)*ratio*mvnpdf(UpdateBETA',MyEvn.mean',MyEvn.var)*InvGDen(UpdateSigma2,MyEvn.A,MyEvn.B);
    %ARatio=(1-Update.d)/(MyEvn.n-Old.Num)*ratio;
    PRatio=(MyEvn.n-Old.Num)/BSDensity(UpdateBETA,UpdateSigma2,MyEvn.X(newIndex,:),MySim.obsNu(newIndex)',Old.sigma2,MyEvn);
    %PRatio=MyEvn.n-Old.Num;
    multi=Like*ARatio*PRatio;
    alpha=min(multi,1);
    u=rand(1);
    if u <= alpha
        New=Update;
    end
end
bloop=bloop+1;
end

function [ New,alpha,bloop,multi ] = ShiftMove( Old,MyEvn,MySim,bloop )
W = MyEvn.W;
Update = Old;
indic = NaN(1,Old.Num);
for i = 1:Old.Num
    Index = find(MyEvn.DistIndex==Old.Center(i));
    group = MyEvn.DistIndex(W(Index,:)==1);
    for j = 1:Old.Num
        group(group==Old.Center(j)) = [];
    end
    gsize = length(group);
    if gsize==0
        indic(i) = 0;
    else indic(i) = 1;
    end
end
subCenter = Old.Center(indic==1);
nsize = length(subCenter);
target = subCenter(randsample(nsize,1));
Index = find(MyEvn.DistIndex==target);
group = MyEvn.DistIndex(W(Index,:)==1);
for i=1:Old.Num
    group(group==Old.Center(i)) = [];
end
gsize=length(group);
if gsize==1
    newCenter=group;
elseif gsize > 1
    newCenter=group(randsample(gsize,1));
else
    disp('empty!')
end
Update.Center(Old.Center==target)=newCenter;
Update.Cluster=ClusterGen(MyEvn.DistMatr,Update.Center,MyEvn.DistIndex);
count = histc(Update.Cluster, sort(unique(Update.Cluster)));
New = Old; alpha = -99; multi = -99;
if min(count) >= MyEvn.nrmin
    ARatio=NuDensity(Update.Nu,Update.Center,Update.Cluster,Update.BETA,Update.Sigma2,Update.sigma2,MyEvn.X)/NuDensity(Old.Nu,Old.Center,Old.Cluster,Old.BETA,Old.Sigma2,Old.sigma2,MyEvn.X);
    indic=NaN(1,Update.Num);
    for i=1:Update.Num
        Index=find(MyEvn.DistIndex==Update.Center(i));
        group=MyEvn.DistIndex(W(Index,:)==1);
        for j=1:Update.Num
            group(group==Update.Center(j))=[];
        end
        gsize=length(group);
        if gsize==0
            indic(i)=0;
        else indic(i)=1;
        end
    end
    if sum(indic)==0
        New=Old;
        alpha=0;
        multi=0;
    else
        subCenter=Update.Center(indic==1);
        Revnsize=length(subCenter);
        target=newCenter;
        Index=find(MyEvn.DistIndex==target);
        group=MyEvn.DistIndex(W(Index,:)==1);
        for i=1:Old.Num
            group(group==Update.Center(i))=[];
        end
        Revgsize=length(group);
        PRatio=nsize*gsize/(Revnsize*Revgsize);
        multi=ARatio*PRatio;
        alpha=min(1,multi);
        u=rand(1);
        if u<=alpha
            New=Update;
        end
    end
end
bloop=bloop+1;
end

function [ NewSigma2,NewBETA ] = UpdateNew( Cluster,Center,OldResult,MyEvn,MySim )
trialnum=length(MyEvn.Sigma2Range);
Sigma2Range=MyEvn.Sigma2Range;
Index=find(Cluster==Center);
nnew=length(Index);
Xnew=MyEvn.X(Index,:);
%Nunew=OldResult.Nu(Index);
Nunew=MySim.obsNu(Index)';
Expectednew=MyEvn.Expected(Index);
Obsnew=MySim.Obs(Index);
loglike=NaN(1,trialnum);
for i=1:trialnum
    SIGMAinv=1/Sigma2Range(i)*eye(nnew)-OldResult.sigma2/Sigma2Range(i)/(Sigma2Range(i)+nnew*OldResult.sigma2)*ones(nnew,nnew);
    A=Xnew'*SIGMAinv*Xnew+MyEvn.invvar;
    B=Xnew'*SIGMAinv*Nunew'+MyEvn.invvar*MyEvn.mean;
    R=Nunew*SIGMAinv*Nunew'+MyEvn.mean'*MyEvn.invvar*MyEvn.mean;
    loglike(i)=-(MyEvn.A+(nnew+1)/2)*log(Sigma2Range(i))-MyEvn.B/Sigma2Range(i)-0.5*log(Sigma2Range(i)+nnew*OldResult.sigma2)-0.5*log(det(A))-0.5*R+0.5*B'*inv(A)*B;
end
Maxloglike=max(loglike);
w=exp(loglike-Maxloglike)/sum(exp(loglike-Maxloglike));
u = rand(1);
cump = cumsum([0 w(1:(end-1))]);
i0 = sum(u > cump);
Sigma2upd=Sigma2Range(i0);
SIGMAinv=1/Sigma2upd*eye(nnew)-OldResult.sigma2/Sigma2upd/(Sigma2upd+nnew*OldResult.sigma2)*ones(nnew,nnew);
A=Xnew'*SIGMAinv*Xnew+MyEvn.invvar;
B=Xnew'*SIGMAinv*Nunew'+MyEvn.invvar*MyEvn.mean;
Ainv=inv(A);
nsize=size(A,2);
for m=1:(nsize-1)
    for n=(m+1):nsize
        Ainv(m,n)=Ainv(n,m);
    end
end
BETAupd=mvnrnd(Ainv*B,Ainv);
Qinv=inv(SIGMAinv+diag(Expectednew'.*exp(Nunew)));
for m=1:(nnew-1)
    for n=(m+1):nnew
        Qinv(m,n)=Qinv(n,m);
    end
end
NewSigma2=Sigma2upd;
NewBETA=BETAupd';
end

function [New,alpha,bloop,multi,ARatio,PRatio ] = DeathMove( Old,MyEvn,MySim,bloop )
DelCenter=randsample(Old.Center,1);
Update.Center=Old.Center;
Update.Center(Old.Center==DelCenter)=[];
Update.Cluster=ClusterGen(MyEvn.DistMatr,Update.Center,MyEvn.DistIndex);
count = histc(Update.Cluster, sort(unique(Update.Cluster)));
New = Old; alpha = -99; multi = -99; ARatio = -99; PRatio = -99;
if min(count) >= MyEvn.nrmin
    Update.Num=length(Update.Center);
    DelIndex=find(Old.Center==DelCenter);
    Update.BETA=Old.BETA(:,setdiff(1:end,DelIndex));
    Update.Sigma2=Old.Sigma2(setdiff(1:end,DelIndex));
    Update.sigma2=Old.sigma2;
    Update.d=Old.d;
    Update.Nu=Old.Nu;
    UpdateIndex=find(Old.Cluster==DelCenter);
    %[Update.Nu UpdateProp]=DeleteOld(Update,DelCenter,Old,MyEvn,MySim);
    Like=1;
    ratio=NuDensity(Update.Nu,Update.Center,Update.Cluster,Update.BETA,Update.Sigma2,Update.sigma2,MyEvn.X)/NuDensity(Old.Nu,Old.Center,Old.Cluster,Old.BETA,Old.Sigma2,Old.sigma2,MyEvn.X);
    ARatio=(MyEvn.n-Update.Num)/(1-Update.d)*ratio/mvnpdf(Old.BETA(:,DelIndex)',MyEvn.mean',MyEvn.var)/InvGDen(Old.Sigma2(DelIndex),MyEvn.A,MyEvn.B);
    %ARatio=(MyEvn.n-Update.Num)/(1-Update.d)*ratio;
    %ARatio=(MyEvn.n-Update.Num)/(1/(Update.d-(1-sqrt(5))/2)+(1-sqrt(5))/2)*ratio/mvnpdf(Old.BETA(:,DelIndex)',MyEvn.mean',MyEvn.var)/InvGDen(Old.Sigma2(DelIndex),MyEvn.A,MyEvn.B);
    PRatio=1/(MyEvn.n-Update.Num)*BSDensity(Old.BETA(:,DelIndex),Old.Sigma2(DelIndex),MyEvn.X(UpdateIndex,:),MySim.obsNu(UpdateIndex)',Update.sigma2,MyEvn);
    %PRatio=1/(MyEvn.n-Update.Num);
    multi=Like*ARatio*PRatio;
    alpha=min(multi,1);
    u=rand(1);
    if u<=alpha
        New=Update;
    end
end
bloop=bloop+1;
end

function [ New,alpha,bloop,multi ] = SwitchMove( Old,MyEvn,MySim,bloop )
index=randsample(Old.Num,2);
Update=Old;
Update.Center(index(1))=Old.Center(index(2));
Update.Center(index(2))=Old.Center(index(1));
Update.Cluster=ClusterGen(MyEvn.DistMatr,Update.Center,MyEvn.DistIndex);
Update.BETA(:,index(1))=Old.BETA(:,index(2));
Update.BETA(:,index(2))=Old.BETA(:,index(1));
Update.Sigma2(index(1))=Old.Sigma2(index(2));
Update.Sigma2(index(2))=Old.Sigma2(index(1));
ARatio=NuDensity(Update.Nu,Update.Center,Update.Cluster,Update.BETA,Update.Sigma2,Update.sigma2,MyEvn.X)/NuDensity(Old.Nu,Old.Center,Old.Cluster,Old.BETA,Old.Sigma2,Old.sigma2,MyEvn.X);
multi=ARatio;
alpha=min(1,multi);
u=rand(1);
if u<=alpha
    New=Update;
else New=Old;
end
bloop=bloop+1;
end

function [ New,bloop ] = HyperMove( Old,MyEvn,MySim,bloop )
New=Old;
sigma2=Old.sigma2;
invvar=MyEvn.invvar;
range=MyEvn.Sigma2Range;
trial=length(range);
for i=1:Old.Num
    Sigma2=Old.Sigma2(i);
    Index=find(Old.Cluster==Old.Center(i));
    X=MyEvn.X(Index,:);
    Nu=Old.Nu(Index);
    nsub=length(Index);
    SIGMAinv=1/Sigma2*eye(nsub)-sigma2/Sigma2/(Sigma2+nsub*sigma2)*ones(nsub,nsub);
    A=X'*SIGMAinv*X+invvar;
    B=X'*SIGMAinv*Nu'+invvar*MyEvn.mean;
    invA=inv(A);
    nsize=size(A,2);
    for m=1:(nsize-1)
        for n=(m+1):nsize
            invA(m,n)=invA(n,m);
        end
    end
    New.BETA(:,i)=mvnrnd(invA*B,invA)';
    BETA=New.BETA(:,i);
    loglike=NaN(1,trial);
    for j=1:trial
        SIGMAinv=1/range(j)*eye(nsub)-sigma2/range(j)/(range(j)+nsub*sigma2)*ones(nsub,nsub);
        loglike(j)=-(MyEvn.A+(nsub+1)/2)*log(range(j))-MyEvn.B/range(j)-0.5*log(range(j)+nsub*sigma2)-0.5*((Nu'-X*BETA)'*SIGMAinv*(Nu'-X*BETA));
    end
    Maxloglike=max(loglike);
    w=exp(loglike-Maxloglike)/sum(exp(loglike-Maxloglike));
    u = rand(1);
    cump = cumsum([0 w(1:(end-1))]);
    i0 = sum(u > cump);
    New.Sigma2(i)=range(i0);
end
loglike=NaN(1,trial);
for i=1:trial
    loglike(i)=-(MyEvn.a+1)*log(range(i))-MyEvn.b/range(i);
    for j=1:New.Num
        Sigma2=New.Sigma2(j);
        Index=find(New.Cluster==New.Center(j));
        nsub=length(Index);
        X=MyEvn.X(Index,:);
        Nu=New.Nu(Index);
        BETA=New.BETA(:,j);
        loglike(i)=loglike(i)-0.5*log(Sigma2+nsub*range(i))+range(i)/2/Sigma2/(Sigma2+nsub*range(i))*(Nu'-X*BETA)'*ones(nsub,nsub)*(Nu'-X*BETA);
    end
end
Maxloglike=max(loglike);
w=exp(loglike-Maxloglike)/sum(exp(loglike-Maxloglike));
u = rand(1);
cump = cumsum([0 w(1:(end-1))]);
i0 = sum(u > cump);
New.sigma2=range(i0);
bloop=bloop+1;
end

function [ New,alpha,bloop,multi ] = HeightMove( Old,MyEvn,MySim,bloop )
New=Old;
multiple=NaN(1,Old.Num);
Alpha=NaN(1,Old.Num);
sigma2=Old.sigma2;
for i=1:Old.Num;
    Sigma2=Old.Sigma2(i);
    BETA=Old.BETA(:,i);
    Index=find(Old.Cluster==Old.Center(i));
    nsub=length(Index);
    X=MyEvn.X(Index,:);
    Exp=MyEvn.Expected(Index);
    Obs=MySim.Obs(Index);
    SIGMAinv=1/Sigma2*eye(nsub)-sigma2/Sigma2/(Sigma2+nsub*sigma2)*ones(nsub,nsub);
    Qinv=inv(SIGMAinv+diag(Exp'.*exp(Old.Nu(Index))));
    for m=1:(nsub-1)
        for n=(m+1):nsub
            Qinv(m,n)=Qinv(n,m);
        end
    end
    P=Obs'-Exp.*exp(Old.Nu(Index)')+Exp.*(Old.Nu(Index)').*exp(Old.Nu(Index)');
    MU=Qinv*(SIGMAinv*X*BETA+P);
    Nuupd=mvnrnd(MU',Qinv);
    UpdateProp=mvnpdf(Nuupd,MU',Qinv);
    Like=prod(exp((Nuupd-Old.Nu(Index)).*Obs-Exp'.*(exp(Nuupd)-exp(Old.Nu(Index)))));
    mu=X*BETA;
    SIGMA=Sigma2*eye(nsub)+sigma2*ones(nsub,nsub);
    ARatio=mvnpdf(Nuupd,mu',SIGMA)/mvnpdf(Old.Nu(Index),mu',SIGMA);
    PRev=Obs'-Exp.*exp(Nuupd')+Exp.*(Nuupd').*exp(Nuupd');
    MURev=Qinv*(SIGMAinv*X*BETA+PRev);
    RevProp=mvnpdf(Old.Nu(Index),MURev',Qinv);
    PRatio=RevProp/UpdateProp;
    multiple(i)=Like*ARatio*PRatio;
    Alpha(i)=min(multiple(i),1);
    u=rand(1);
    if u<=Alpha(i)
        New.Nu(Index)=Nuupd;
    end
end
bloop=bloop+1;
multi=mean(multiple);
alpha=mean(Alpha);
end

function [ p ] = NuDensity( Nu,Center,Cluster,BETA,Sigma2,sigma2,X )
n=length(Center);
p=1;
for i=1:n
    Index=find(Cluster==Center(i));
    ncluster=length(Index);
    MU=X(Index,:)*BETA(:,i);
    SIGMA=Sigma2(i)*eye(ncluster)+sigma2*ones(ncluster,ncluster);
    p=p*mvnpdf(Nu(Index),MU',SIGMA);
end
end

function [ DW ] = DistMatrix(Latitude,Longitude)
n=length(Latitude);
DW=NaN(n,n);
for i=1:n
    for j=1:n
        DW(i,j)=sqrt((Latitude(i)-Latitude(j))^2+(Longitude(i)-Longitude(j))^2);
    end
end
end

function [ Cluster ] = ClusterGen( DistMatr,CenterIndex,DistIndex )
n=length(DistIndex);
Cluster=DistIndex;
Index=NaN(1,length(CenterIndex));
for i=1:length(CenterIndex)
    Index(i)=find(DistIndex==CenterIndex(i));
end
for i=setdiff(1:n,Index)
    %Cluster(i)=DistMin(Coordinate(i,:),Coordinate(Index,:),CenterIndex);
    Cluster(i)=DistIndex(DistMin2(DistMatr,i,Index));
end
end

function [ MinIndex ] = DistMin2( DistMatr,DistIndex,CenterIndex)
%Here DistIndex and CenterIndex must be the order index
DistCand=DistMatr(DistIndex,CenterIndex);
MinIndex=find(DistCand==min(DistCand));
MinIndex=CenterIndex(MinIndex(1));
end

function [ p ] = InvGDen( x,a,b )
p=b^a/gamma(a)*x^(-a-1)*exp(-b/x);
end

function [ prob ] = BSDensity( beta,Sigma2,X,Nu,sigma2,MyEvn )
p=size(MyEvn.X,2);
n=length(Nu);
par=0.001:0.001:1;
trial=length(par);
SIGMAinv=1/Sigma2*eye(n)-sigma2/Sigma2/(Sigma2+n*sigma2)*ones(n,n);
invvar=inv(MyEvn.var);
A=X'*SIGMAinv*X+invvar;
B=X'*SIGMAinv*Nu'+invvar*MyEvn.mean;
R=Nu*SIGMAinv*Nu'+MyEvn.mean'*invvar*MyEvn.mean;
Norm=0;
for i=1:(trial-1)
    Norm=Norm+0.001/6*(SDensity(par(i),sigma2,MyEvn,A,B,R,p,n)+4*SDensity((par(i)+par(i+1))/2,sigma2,MyEvn,A,B,R,p,n)+SDensity(par(i+1),sigma2,MyEvn,A,B,R,p,n));
end
prob=1/Norm/(Sigma2^((n-1)/2))/(Sigma2+n*sigma2)^0.5*exp(-0.5*(Nu'-X*beta)'*SIGMAinv*(Nu'-X*beta))*exp(-0.5*(beta-MyEvn.mean)'*inv(MyEvn.var)*(beta-MyEvn.mean))/(Sigma2^(MyEvn.A+1))*exp(-MyEvn.B/Sigma2);
end

function [ prob ] = SDensity( Sigma2,sigma2,MyEvn,A,B,R,p,n )
prob=(sqrt(2*pi))^p/(Sigma2^(MyEvn.A+(n+1)/2))*exp(-MyEvn.B/(Sigma2))/((Sigma2+n*sigma2)^0.5)/(det(A))^0.5*exp(-0.5*R+0.5*B'*inv(A)*B);
end

function [ New,bloop ] = updateD( Old, MyEvn, MySim,bloop )
New=Old;
u=rand(1);
New.d=1-(1-u)^(1/(Old.Num+1));
bloop=bloop+1;
end











clc;clear;
tic;
%数据读入
load('T1.mat');
%%
V1=[100/3.6,10,10,10.02,10.94,14.27];
V2=[100/3.6,9.61,9.39,8.78,7.6,6.92];

mark=["-","-","--",":","--","--"];
Maks=["none",".","none","none","none","none"];
colo=[[0,0,0];[227,140,122];[255,0,0];[0,0,255];[176,101,89];[215,157,164]]./255;
n = length(V1);
for i=1:1:n
    v1=V1(i);
    v2=V2(i);
    %牵引时间计算
    dt=0.01;
    A=zeros(1,1800);
    S1=zeros(1,1800);
    E1=zeros(1,1800);
    for i1=2:1:1800
        f1=2.0895+0.0098*A(i1-1)+0.006*(A(i1-1)^2);
        A(i1)=(Nmax-f1)*dt/(m*rou)+A(i1-1);
        S1(i1)=S1(i1-1)+(A(i1)+A(i1-1))*dt/2;
        E1(i1)=Nmax*A(i1)*dt;
        if(A(i1)>=v1) 
            break;
        end
    end
    A = A(A~=0);
    A=[0,A];
    S1 = S1(S1~=0);
    S1=[0,S1];
    E1 = E1(E1~=0);
    E1=[0,E1];
    t1=(i1-1)*dt;
   %制动时间计算
    D=zeros(1,1800);
    D(1)=v2;
    S4=zeros(1,1800);
    for i4=2:1:1800
        f4=2.0895+0.0098*D(i4-1)+0.006*(D(i4-1)^2);
        D(i4)=-(Stopmax+f4)*dt/(m*rou)+D(i4-1);
        S4(i4)=S4(i4)+(D(i4)+D(i4-1))*dt/2;
        if(D(i4)<=0) 
            break;
        end
    end
    D = D(D~=0);
    S4 = S4(S4~=0);
    S4=[0,S4];
    t4=(i4-1)*dt;
    %惰行
    dt=0.01;
    C=zeros(1,180000);
    C(1)=v1;
    S3=zeros(1,180000);
    for i3=2:1:180000
        f3=2.0895+0.0098*C(i3-1)+0.006*(C(i3-1)^2);
        C(i3)=(-f3)*dt/(m*rou)+C(i3-1);
        S3(i3)=S3(i3-1)+(C(i3)+C(i3-1))*dt/2;
        if(C(i3)<=v2) 
            break;
        end
    end
    C = C(C~=0);
    S3 = S3(S3~=0);
    S3=[0,S3];
    t3=(i3-1)*dt;
    % %计算巡航时间
    S2=S-S1(end)-S4(end)-S3(end);
    t2=S2/v1;
    f2=2.0895+0.0098*v1+0.006*(v1^2);
    S3=S3+S1(end)+S2;
    S4=S4+S3(end);
    E2=f2*v1*dt;
    %速度/距离
    hold on;
    col=colo(i,:);
    mak=mark(i);
    mak2=Maks(i);
    %牵引制动力/距离
%     plot([0,S1(end)],[Nmax,Nmax],'Color',col, 'LineWidth',2);
%     plot([S1(end),S1(end)],[Nmax,f2],'Color',col, 'LineWidth',2);
%     plot([S1(end),S3(1)],[f2,f2],'Color',col, 'LineWidth',2);
%牵引力
    plot([0,S1(end),S1(end),S3(1)],[Nmax,Nmax,f2,f2],'Color',col,  'LineStyle',mak,'Marker',mak2,'LineWidth',2);
%制动力
%     plot([S4(1),S4(end)],[Stopmax,Stopmax],'Color',col,  'LineStyle',mak,'Marker',mak2,'LineWidth',2)
end
legend('200s','210s','220s','250s','350s','500s');
% xlim([xlim, 5144.7]);
% ylim([0,310]);
tt2=toc;
disp(['代码执行时间：', num2str(tt2), '秒']);
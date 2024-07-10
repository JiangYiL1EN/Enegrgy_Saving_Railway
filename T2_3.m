%时间距离
clc;clear;
tic;
%数据读入
load('T1.mat');
%%
EX=[];
V_1X=[];
V_2X=[];
tmin1X=[];
g=9.8;
v1=[100/3.6,23,22,19,14,10];
v2=[86/3.6,23,22,19,14,10];
% v1=100/3.6;
% v2=86/3.6;
v1=17;
v2=15;
mark=["-","-","--",":","--","--"];
Maks=["none",".","none","none","none","none"];
colo=[[0,0,0];[227,140,122];[255,0,0];[0,0,255];[176,101,89];[215,157,164]]./255;
n = length(v1);
for i=1:1:n
    V1=v1(i);
    V2=v2(i);
    
    Emin=80000;
    tmin1=100;

    dt=0.01;
    A=zeros(1,80000);
    S1=zeros(1,80000);
    E1=zeros(1,80000);
    jiao=atan(0.0617284*(10^-3));
    for i1=2:1:80000
        f1=2.0895+0.0098*A(i1-1)+0.006*(A(i1-1)^2);
        %计算ftra
        if  A(i1-1)<=10
            ftra=310;
        else
            ftra=310*(10)/A(i1-1);
        end

        if jiao>=0
            A(i1)=(ftra-(f1+m*g*sin(jiao)))*dt/(m*rou)+A(i1-1);
        else
            A(i1)=(m*g*sin(jiao)+ftra-f1)*dt/(m*rou)+A(i1-1);
        end
        S1(i1)=S1(i1-1)+(A(i1)+A(i1-1))*dt/2;
        E1(i1)=E1(i1-1)+ftra*A(i1)*dt;
        if(A(i1)>V1)
            A(i1)=V1;
        end

        if(S1(i1)>4259.1)&(A(i1)>86/3.6)
            flag=1;
            break;
        end
        if(A(i1)>=V1) 
            break;
        end

        index = find(DisGra(:,1) <S1(i1));
        index=index(end);
        jiao=atan(DisGra(index, 2)*(10^-3));
    end
    if flag==1
        continue;
    end
    A = A(A~=0);
    A=[0,A];
    S1 = S1(S1~=0);
    S1=[0,S1];
    E1 = E1(E1~=0);
    E1=[0,E1];
    t1=(i1-1)*dt;

    %制动时间计算
    D=zeros(1,80000);
    % D(1)=Vmax;
    S4=zeros(1,80000);
    E4=zeros(1,80000);
    jiao2=atan(20.3086*(10^-3));
    for i4=2:1:80000
        f4=2.0895+0.0098*D(i4-1)+0.006*(D(i4-1)^2);
        if  D(i4-1)<=17
            ftra=260;
        else
            ftra=260*(17)/D(i4-1);
        end

        if jiao2>=0
            D(i4)=-(0-(f4+m*g*sin(jiao2)+ftra))*dt/(m*rou)+D(i4-1);
        else
            D(i4)=-(m*g*sin(jiao2)-(ftra+f4))*dt/(m*rou)+D(i4-1);
        end

        S4(i4)=S4(i4)+(D(i4)+D(i4-1))*dt/2;
        E4(i4)=E4(i4-1)+ftra*D(i4)*dt;
        if(D(i4)>V2) 
            break;
        end

        index = find(DisGra(:,1) <S-S4(i4));
        index=index(end);
        jiao2=atan(DisGra(index, 2)*(10^-3));
    end
    D = D(D~=0);
    D=[0,D];
    D=flip(D);
    S4 = S4(S4~=0);
    S4=[0,S4];
    S4=5144.7-flip(S4);
    t4=length(D)*dt;

    %惰行
    dt=0.01;
    C=zeros(1,180000);
    C(1)=V2;
    S3=zeros(1,180000);
    for i3=2:1:180000
        f3=2.0895+0.0098*C(i3-1)+0.006*(C(i3-1)^2);

        if jiao2>=0
            C(i3)=-(0-(f3+m*g*sin(jiao2)))*dt/(m*rou)+C(i3-1);
        else
            C(i3)=-(m*g*sin(jiao2)-(f3))*dt/(m*rou)+C(i3-1);
        end
       
        S3(i3)=S3(i3-1)+(C(i3)+C(i3-1))*dt/2;
        
        if(S4(1)-S3(i3)>4259.1)&(C(i3)>86/3.6)
            C(i3)=86/3.6;
            break;
        end
        if(C(i3)>=V1) 
            break;
        end
        index = find(DisGra(:,1) <S4(1)-S3(i3));
        index=index(end);
        jiao2=atan(DisGra(index, 2)*(10^-3));
    end
    C = C(C~=0);
    C=flip(C);
    S3 = S3(S3~=0);
    S3=[0,S3];
    S3=S4(1)-flip(S3);
    t3=(i3-1)*dt;

    %巡航
    dt=0.01;
    B=zeros(1,180000);
    B(1)=V1;
    S2=zeros(1,180000);
    S2(1)=S1(end);
    E2=zeros(1,18000);
    E2(1)=E1(end);
    for i2=2:1:180000
        f2=2.0895+0.0098*B(i2-1)+0.006*(B(i2-1)^2);
        B(i2)=B(i2-1);
        if jiao>=0
            ftra=f2+m*g*sin(jiao);
        else
            ftra=f2-m*g*sin(jiao);
        end

        S2(i2)=S2(i2-1)+V1*dt;
        E2(i1)=E2(i1-1)+ftra*B(i1)*dt;
        if(S2(i2)>4259.1)&(B(i2)>86/3.6)
            B(i2)=86/3.6;
            S2(i2)=S2(i2-1)+1;
        end
        if(S2(i2)>(S3(1)))
            break;
        end

        index = find(DisGra(:,1) <S-S2(i2));
        index=index(end);
        jiao=atan(DisGra(index, 2)*(10^-3));
    end

    B = B(B~=0);
    S2 = S2(S2~=0);
    t2=(i2-1)*dt;
    
    hold on;
    col=colo(i,:);
    mak=mark(i);
    mak2=Maks(i);
    plot([S1,S2,S3,S4],[((0:1:i1-1)*dt),(((0:1:i2-1)*dt)+t1),(((0:1:i3-1)*dt)+t1+t2),(((0:1:i4-1)*dt)+t1+t2+t3)],'Color',col, 'LineStyle',mak,'Marker',mak2, 'LineWidth',2)
end
legend('206s','216s','226s','256s','356s','506s','FontSize', 12);
set(gca, 'FontSize', 12);
xlim([0, 5200]);
% ylim([0,30]);

tt2=toc;
disp(['代码执行时间：', num2str(tt2), '秒']);
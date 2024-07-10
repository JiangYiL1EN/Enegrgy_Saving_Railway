%最短时间
clc;clear;
%数据读入
load('T1.mat');
Vtra=10;
Ftra=310;
Vtra2=17;
Ftra2=260;
g=9.8;
%%
%牵引时间计算
dt=0.01;
A=zeros(1,18000);
S1=zeros(1,18000);
E1=zeros(1,18000);
jiao=0.0617284;
for i1=2:1:18000
    f1=2.0895+0.0098*A(i1-1)+0.006*(A(i1-1)^2);
    %计算ftra
    if  A(i1-1)<=Vtra
        ftra=Ftra;
    else
        ftra=Ftra*(Vtra)/A(i1-1);
    end
    
    if jiao>=0;
        A(i1)=(ftra-(f1+m*g*sin(jiao)))*dt/(m*rou)+A(i1-1);
    else
        A(i1)=(m*g*sin(jiao)+ftra-f1)*dt/(m*rou)+A(i1-1);
    end
    S1(i1)=S1(i1-1)+(A(i1)+A(i1-1))*dt/2;
    E1(i1)=E1(i1-1)+ftra*A(i1)*dt;
    
    if(A(i1)>=100/3.6) 
        break;
    end
    
    index = find(DisGra(:,1) <S1(i1));
    index=index(end);
    jiao=DisGra(index, 2);
end
A = A(A~=0);
A=[0,A];
S1 = S1(S1~=0);
S1=[0,S1];
E1 = E1(E1~=0);
E1=[0,E1];
t1=(i1-1)*dt;
%%
%制动时间计算
D=zeros(1,1800);
D(1)=0;
S4=zeros(1,1800);
E4=zeros(1,1800);
jiao2=20.3086;
for i4=2:1:1800
    f4=2.0895+0.0098*D(i4-1)+0.006*(D(i4-1)^2);
    if  D(i4-1)<=Vtra2
        ftra=Ftra2;
    else
        ftra=Ftra2*(Vtra2)/D(i4-1);
    end
    
    if jiao2>=0;
        D(i4)=-(0-(f4+m*g*sin(jiao2)+ftra))*dt/(m*rou)+D(i4-1);
    else
        D(i4)=-(m*g*sin(jiao2)-(ftra+f4))*dt/(m*rou)+D(i4-1);
    end
    
    S4(i4)=S4(i4-1)+(D(i4)+D(i4-1))*dt/2;
    E4(i4)=E4(i4-1)+ftra*D(i4)*dt;
    if(D(i4)>=84/3.6) 
        break;
    end
    
    index = find(DisGra(:,1) <S-S4(i4));
    index=index(end);
    jiao2=DisGra(index, 2);
end
D = D(D~=0);
S4 = S4(S4~=0);
S4=[0,S4];
S4=5144.7-flip(S4);
t4=(i4-1)*dt;
%%
%计算巡航时间
S2=4258.1-S1(end);
t2=(4258.1-S1(end))/Vmax+(S4(1)-4259.1)/(86/3.6);
%%
%整合
T=t1+t2+t4+1
disp(T);
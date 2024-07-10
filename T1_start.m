%计算最短时间
clc;clear;
%数据读入
load('T1.mat');
%%
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
    if(A(i1)>=Vmax) 
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
%%
%制动时间计算
D=zeros(1,1800);
D(1)=Vmax;
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
%%
% %计算巡航时间
S2=S-S1(end)-S4(end);
t2=S2/Vmax;
f2=2.0895+0.0098*Vmax+0.006*(Vmax^2);
S4=S1(end)+S2+S4;
E2=f2*Vmax*dt;
%%
%整合
T=t1+t2+t4;
disp(T);
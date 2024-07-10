clc;clear;
%数据读入
load('T1.mat');
%%
EX=[];
V11X=[];
V21X=[];
tmin1X=[];
for ZZ=[210,220,250,350,500]%
    Emin=800000;
    for V1=20%
        tmin=50;
        for V2=0:0.01:V1
            %牵引时间计算
            dt=0.01;
            A=zeros(1,1800);
            S1=zeros(1,1800);
            E1=zeros(1,1800);
            for i1=2:1:1800
                f1=2.0895+0.0098*A(i1-1)+0.006*(A(i1-1)^2);
                A(i1)=(Nmax-f1)*dt/(m*rou)+A(i1-1);
                S1(i1)=S1(i1-1)+(A(i1)+A(i1-1))*dt/2;
                E1(i1)=E1(i1-1)+Nmax*A(i1)*dt;
                if(A(i1)>=V1) 
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
            D(1)=V2;
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
            C(1)=V1;
            S3=zeros(1,180000);
            for i3=2:1:180000
                f3=2.0895+0.0098*C(i3-1)+0.006*(C(i3-1)^2);
                C(i3)=(-f3)*dt/(m*rou)+C(i3-1);
                S3(i3)=S3(i3-1)+(C(i3)+C(i3-1))*dt/2;
                if(C(i3)<=V2) 
                    break;
                end
            end
            C = C(C~=0);
            S3 = S3(S3~=0);
            S3=[0,S3];
            t3=(i3-1)*dt;

            S2=S-S1(end)-S3(end)-S4(end);
            if S2<0
               continue; 
            end

            t2=S2/Vmax;
            f2=2.0895+0.0098*Vmax+0.006*(Vmax^2);
            S4=S1(end)+S2+S4;
            E2=f2*Vmax*t2;

            if abs(((t1+t2+t3+t4)-ZZ))<abs(tmin)
                tmin=(t1+t2+t3+t4)-ZZ;
                E=E2+E1(end);
                V22=V2;
            end
        end
        if(E<Emin)
            Emin=E;
            V11=V1;
            V21=V22;
            tmin1=tmin;
        end
        
    end
    EX=[EX,Emin];
    V11X=[V11X,V11];
    V21X=[V21X,V21];
    tmin1X=[tmin1X,tmin1];
end
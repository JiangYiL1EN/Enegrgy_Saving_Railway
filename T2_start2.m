clc;clear;
%数据读入
load('T1.mat');
%%
EX=[];
V_1X=[];
V_2X=[];
tmin1X=[];
g=9.8;
for ZZ=[10,20,50,150,300]+207.34%,
    Emin=80,000;
    tmin1=100;
    for V1=1:0.01:100/3.6%
        t_lower=8;
        E=0;
%         V22=3;
        for V2=1:0.01:min(V1,86/3.6)%
            flag=0;
            %牵引时间计算
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
%                 if(S4(i4)>4259.1)&(D(i4)>=86/3.6)
%                     flag=1;
%                     break;
%                 end
                if(D(i4)>V2) 
                    break;
                end

                index = find(DisGra(:,1) <S-S4(i4));
                index=index(end);
                jiao2=atan(DisGra(index, 2)*(10^-3));
            end
%             if flag==1
%                 continue;
%             end
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
                
% %                 if(S3(i3)>4259.1)&(C(i3)>=86/3.6)
% %                     flag=1;
% %                     break;
% %                 end

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
%             if flag==1
%                 continue;
%             end
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
                if(S2(i2)>(S3(1)-S1(end)))
                    break;
                end

                index = find(DisGra(:,1) <S-S2(i2));
                index=index(end);
                jiao=atan(DisGra(index, 2)*(10^-3));
            end
            if flag==1
                continue;
            end
            B = B(B~=0);
            S2 = S2(S2~=0);
            t2=(i2-1)*dt;
            
            %判断
            if abs(((t1+t2+t3+t4)-ZZ))<abs(t_lower)
                t_lower=(t1+t2+t3+t4)-ZZ;
                E=E2(end)+E1(end)+E4(end);
                V22=V2;
            end
        end
        if(E<Emin)&(E~=0)
            Emin=E;
            V_1=V1;
            V_2=V22;
            tmin1=t_lower;
        end 
%        if abs(t_lower)<abs(tmin1)
%             Emin=E;
%             V_1=V1;
%             V_2=V22;
%             tmin1=t_lower;
%        end
    end
    EX=[EX,Emin];
    V_1X=[V_1X,V_1];
    V_2X=[V_2X,V_2];
    tmin1X=[tmin1X,tmin1];
end
disp(V_1X);
disp(V_2X);
clc
close all
clear all
nu=3*10^(-6);
U=30;
eta=(0:0.1:10);
x=(0:0.1:10);
N=length(eta);
F=zeros(N,5);
h=0.1;
pr=1;
er=0.001; 
IC1=0;
ic=[0 0 IC1];
f1=rk4sys(eta,ic);
IC2=1;
ic=[0 0 IC2];
f2=rk4sys(eta,ic);
%err=er+1;
f=@(eta,y) [y(2) y(3) -y(1)*y(3)  y(5) -pr*y(1)*y(5)];
    err=1;
    err2=1;
    tol=1e-4;
    bv=[0.1 0.3];
    bv2=bv;
        while err>tol || err2>tol
            for j=1:2
                F=zeros(N,5);
                F(1,:)=[0;0;bv(j);0;bv2(j)];

                for i = 1:N-1
                    K1 = h * f(eta(i) , F(i,:)) ; 
                    K2 = h * f(eta(i)+h/2 , F(i,:)+K1/2) ; 
                    K3 = h * f(eta(i)+h/2 , F(i,:)+K2/2) ; 
                    K4 = h * f(eta(i)+h , F(i,:)+K3) ; 
                    F(i+1,:) = F(i,:) + 1/6*(K1 + 2*K2 + 2*K3+K4) ;
                end
                g_inf(j) = F(end,2);
                g_inf2(j) = F(end,4);
            end
        [err,index]= max(abs(1-g_inf));
        [err2,index2]= max(abs(1-g_inf2));

        P = diff(g_inf);
        P2 = diff(g_inf2);
    
    

        if P~0;
            bv_new = bv(1)+(diff(bv)/diff(g_inf))*(1-g_inf(1)); %interpolated value
            bv(index) = bv_new; 
        end
        if P2~0;
            bv2_new = bv2(1)+(diff(bv2)/diff(g_inf2))*(1-g_inf2(1)); %interpolated value
            bv2(index2) = bv2_new ;
        end
        end
        %{
while err>er
      IC3=((IC2-IC1)/(f2(end,2)-f1(end,2)))*(1-f1(end,2))+IC1;
      ic=[0 0 IC3];
      f3=rk4sys(eta,ic);
      err=abs(f3(end,2)-1);
      IC1=IC2;
      IC2=IC3;
      f1=f2;
      f2=f3;
end
        %}
f=F(:,1);
f_prime=F(:,2);
subplot(1,3,1)
plot(f_prime,eta,'b','linewidth',2)
grid
xlabel('u = df / d\eta')
ylabel('\eta')
%n=length(eta);
%del_eta=eta(2)-eta(1);
%{
for Pr=1
    A=zeros(n,n);
    B=zeros(n,1);
    A(1,1)=1;   
    A(n,n)=1;
    B(n)=1;
    for i=2:n-1
        A(i,i)=-2;
        A(i,i-1)=1-(1/2)*Pr*f(i)*del_eta;
        A(i,i+1)=1+(1/2)*Pr*f(i)*del_eta;
    end
    theta=A\B;
    subplot(1,3,2)
    hold on
    plot(theta,eta,'r','linewidth',1.5)
end 
%}
theta=F(:,4);
subplot(1,3,2)
hold on
plot(theta,eta,'r','linewidth',1.5)
grid
xlabel('\theta')
ylabel('\eta')
legend('Pr = 1','Location','northwest')
u=f_prime*U;

for i=1:length(x)
    for j=1:length(eta)
    Y(i,j) =eta(j).*sqrt(2*nu/U).*sqrt(x(i));
    end
end
for i=1:length(x)
    for j=1:length(eta)
    psi(i,j) =f(j).*sqrt(2*nu*U).*sqrt(x(i));
    end
end
for i=1:length(x)
    for j=1:length(eta)
    u_2D(i,j) =f_prime(j)*U;
    end
end
X=ones(length(x)).*x';

v=sqrt(nu*U./(2.*x)).*(eta.*f_prime'-f');
T=theta*(-35)+60;
for i=1:length(x)
    for j=1:length(eta)
    v_2D(i,j) =sqrt(nu*U./(2.*x(i))).*(eta(j).*f_prime(j)-f(j));
    end
end


figure(2)
   contour(X,Y,psi,50)
   hold on
   plot(X(:,36), Y(:,36));
   
   
  
   
   %y_bl=sqrt(2.*nu.*x./U).*3.5;
   y_bl=linspace(0,0.005,101);
   [XX YY]=meshgrid(x,y_bl);
   
   for i=1:length(x)
       for j=1:length(y_bl)
           eta_new(i,j)=y_bl(j)*sqrt(U/(2*nu*x(i)));
           
       end
   end
  theta_prime=F(:,5);
  theta_prime(1)=0.332*pr^(1/3);
   for i=1:size(eta_new,1)
       for j=1:size(eta_new,2)
           n=round((eta_new(i,j)-mod(eta_new(i,j),0.1))/0.1);
           nn(i,j)=n;
           
           
           if i==1
               u_mesh(i,j)=0;
               v_mesh(i,j)=0;
               T_mesh(i,j)=25;
           elseif j==1
               u_mesh(i,j)=U;
               v_mesh(i,j)=0;
               T_mesh(i,j)=60;
           elseif n>=size(eta_new,1)
               
               u_mesh(i,j)=u(end);
               v_mesh(i,j)=v(end);
               T_mesh(i,j)=T(end);
               f_mesh(i,j)=f(end);
               theta_mesh(i,j)=theta(end);
               theta_prime_mesh(i,j)=theta_prime(end);
           elseif n>1 && n<=size(eta_new,1)-1
               u_mesh(i,j)=u(n+1)+ (u(n+1)-u(n))/(0.1)*(eta_new(i,j)-eta(n));
               v_mesh(i,j)=v(n+1)+ (v(n+1)-v(n))/(0.1)*(eta_new(i,j)-eta(n));
               T_mesh(i,j)=T(n+1)+ (T(n+1)-T(n))/(0.1)*(eta_new(i,j)-eta(n));
               f_mesh(i,j)=f(n+1)+ (f(n+1)-f(n))/(0.1)*(eta_new(i,j)-eta(n));
               theta_mesh(i,j)=theta(n+1)+ (theta(n+1)-theta(n))/(0.1)*(eta_new(i,j)-eta(n));
               theta_prime_mesh(i,j)=theta_prime(n+1)+ (theta_prime(n+1)-theta_prime(n))/(0.1)*(eta_new(i,j)-eta(n));
               
           end
       end
   end
   for i=1:size(eta_new,1)
       for j=1:size(eta_new,2)
            n=round((eta_new(i,j)-mod(eta_new(i,j),0.1))/0.1);
            H(i,j)=-x(i)^0.5*(f_mesh(i,j)*(theta_mesh(i,j)-1)+2*theta_prime_mesh(i,j)/pr);
       end
   end
   
   
  H(:,1)=-2./pr.*theta_prime(1).*x.^0.5;
  
  figure(3)
  contour(XX,YY,H,50)
  
   
function Y=rk4sys(X,ic)
    n=length(X);
    a=X(1);
    b=X(end);
    h=(b-a)/(n-1);
    Y(1,:)=ic';
    y=ic;
    for i=1:n-1
        K1=h*rk_fun(X(i),y);
        K2=h*rk_fun(X(i)+h/2,y+K1/2);
        K3=h*rk_fun(X(i)+h/2,y+K2/2);
        K4=h*rk_fun(X(i)+h,y+K3);
        y=y+(K1+2*K2+2*K3+K4)/6;
        Y(i+1,:)=y';
    end
end
function du=rk_fun(x,u)
    du(1)=u(2);
    du(2)=u(3);
    du(3)=-u(1)*u(3);
end


    
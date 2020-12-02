clc
clear all
close all

 font = 'Times New Roman';
%% run Matrix_H_6pmu.m 
% state estimation and spoofing theta detection

Nb=14;
Npmu=6;
a=[0;1 ;0 ;1; 0; 1; 1; 0; 0; 1; 0; 0;0 ; 1] ; %

%% lines connected to bus n
n1=4  ;n2=5 ;n3=4 ; n4=3 ;  n5=2 ;  n6=2  ;
Mn=[2*(n1+1),2*(n2+1),2*(n3+1),2*(n4+1),2*(n5+1),2*(n6+1)];
%% number of measurements

for i=1:1:Npmu
    sm=sum(Mn(1:i-1));
    num_meas=sum(Mn);
end
%% sigma
sg1=0.02 ;  sg2=0.03 ;sg3=0.03 ;sg4=0.02 ;  sg5=0.01 ; sg6=0.02 ;  
Sg=[sg1,sg2,sg3,sg4,sg5,sg6];
sig=zeros(num_meas,1);

for k=1:1:Npmu
    sm=sum(Mn(1:k-1));
    sig(sm+1:sm+Mn(k))=Sg(k)*ones(Mn(k),1); 
end
Sigma=diag(sig);
invers_sigm = inv(Sigma);
%% W
mu=0;
W=zeros(num_meas,1);
for i=1:1:Npmu
    sm=sum(Mn(1:i-1));
    W(sm+1:sm+Mn(i))=normrnd(mu,Sg(i),Mn(i),1);
end

%% H
H_matrix_6pmu
S_states_14bus
S;
 
V=[real(S(1,1));imag(S(1,1));real(S(2,1));imag(S(2,1));real(S(3,1));imag(S(3,1));real(S(4,1));imag(S(4,1));real(S(5,1));imag(S(5,1));real(S(6,1));imag(S(6,1));real(S(7,1));imag(S(7,1))
    real(S(8,1));imag(S(8,1));real(S(9,1));imag(S(9,1));real(S(10,1));imag(S(10,1));real(S(11,1));imag(S(11,1));real(S(12,1));imag(S(12,1));real(S(13,1));imag(S(13,1));real(S(14,1));imag(S(14,1))];
%% %% phasor measurements
Z_true=H*V;

%% Gama (true)  attack is occured on pmus

delta_theta_spf1=input('delta_teta_spf1=')
delta_theta_spf2=input('delta_teta_spf2=')
 delta_theta_spf3=input('delta_teta_spf3=')
delta_theta_spf4=input('delta_teta_spf4=')
delta_theta_spf5=input('delta_teta_spf5=')
delta_theta_spf6=input('delta_teta_spf6=')
gama1=[cos(delta_theta_spf1),-sin(delta_theta_spf1);sin(delta_theta_spf1),cos(delta_theta_spf1)];
gama2=[cos(delta_theta_spf2),-sin(delta_theta_spf2);sin(delta_theta_spf2),cos(delta_theta_spf2)];
gama3=[cos(delta_theta_spf3),-sin(delta_theta_spf3);sin(delta_theta_spf3),cos(delta_theta_spf3)];
gama4=[cos(delta_theta_spf4),-sin(delta_theta_spf4);sin(delta_theta_spf4),cos(delta_theta_spf4)];
gama5=[cos(delta_theta_spf5),-sin(delta_theta_spf5);sin(delta_theta_spf5),cos(delta_theta_spf5)];
gama6=[cos(delta_theta_spf6),-sin(delta_theta_spf6);sin(delta_theta_spf6),cos(delta_theta_spf6)];
gama=[gama1;gama2;gama3;gama4;gama5;gama6];

Gama_true=eye(num_meas,num_meas);
k=1;
for i=1:1:Nb
  
    if a(i)==1
         
      sm=sum(Mn(1:k-1));

     o=Mn(k);
     j=1;
    
              if i==2
                        while j<=o
              
                          Gama_true(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama1);
                       j=j+2;
                      
                        end
                        k=k+1;
              end
              
              
              
              
              if i==4
                         while j<=o
              
                          Gama_true(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama2);
                          j=j+2;
                         
                         end
                         k=k+1;
              end
              
              if i==6
                         while j<=o
              
                          Gama_true(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama3);
                          j=j+2;
                         
                         end
                         k=k+1;
              end
              
              if i==7
                         while j<=o
              
                          Gama_true(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama4);
                          j=j+2;
                         
                         end
                         k=k+1;
              end
              
               if i==10
                         while j<=o
              
                          Gama_true(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama5);
                          j=j+2;
                         
                         end
                         k=k+1;
               end
              
                if i==14
                         while j<=o
              
                          Gama_true(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama6);
                          j=j+2;
                         
                         end
                         k=k+1;
                end
              
    
  
    
    end
    
end

Gama_true;

Z_at=Gama_true*Z_true+W;

tic
%% gama ( estimation)
%% initialization 
%  theta=0
gam11=1 ;  % cos(theta)
gam12=0 ;  % sin(theta)
gama_new1=[gam11 -gam12;gam12 gam11];
gama_new2=[gam11 -gam12;gam12 gam11];
gama_new3=[gam11 -gam12;gam12 gam11];
gama_new4=[gam11 -gam12;gam12 gam11];
gama_new5=[gam11 -gam12;gam12 gam11];
gama_new6=[gam11 -gam12;gam12 gam11];
Gama_AM=zeros(num_meas,num_meas);
 k=1;
for i=1:1:Nb
     if a(i)==1
         
      sm=sum(Mn(1:k-1));

     o=Mn(k);
     j=1;
              if i==2
                        while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new1);
                       j=j+2;
                      
                        end
                        k=k+1;
              end
              
              
              
              
              if i==4
                         while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new2);
                          j=j+2;
                         
                         end
                         k=k+1;
              end
              
              if i==6
                         while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new3);
                          j=j+2;
                         
                         end
                         k=k+1;
              end
              
              if i==7
                         while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new4);
                          j=j+2;
                         
                         end
                         k=k+1;
              end
              
                if i==10
                         while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new5);
                          j=j+2;
                         
                         end
                         k=k+1;
                end
              
                  if i==14
                         while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new6);
                          j=j+2;
                         
                         end
                         k=k+1;
                  end
             
     end
    
end 

%%  Estimation of 'delta_theta' cuased by spoofing and 'V'

Tolerance= (10)^(-5);
obj=zeros(500,Npmu);
OBJ=zeros(500,1);
count=0;


for d=2:1:500
    
%% step 1  V_hat
G1=[zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb)];

k=1;
for i=1:1:Nb
    if a(i)==1
        sm=sum(Mn(1:k-1));
        G1(:,(k-1)*2*Nb+1:(k-1)*2*Nb+2*Nb)=a(i)*(H(sm+1:sm+Mn(k),:))'*(Gama_AM(sm+1:sm+Mn(k),sm+1:sm+Mn(k)))'* invers_sigm(sm+1:sm+Mn(k),sm+1:sm+Mn(k))*(Gama_AM(sm+1:sm+Mn(k),sm+1:sm+Mn(k)))*H(sm+1:sm+Mn(k),:);
        k=k+1;
    end
end
i=1;
G_AM=zeros(2*Nb,2*Nb);
while i<=Npmu
  
    G_AM=G_AM+G1(:,(i-1)*2*Nb+1:(i-1)*2*Nb+2*Nb);
     i=i+1;
end

k=1;
 SuM1=zeros(2*Nb,2);
for i=1:1:Nb
    if a(i)==1
      sm=sum(Mn(1:k-1));
      SuM1(:,k)=a(i)*(H(sm+1:sm+Mn(k),:))'*(Gama_AM(sm+1:sm+Mn(k),sm+1:sm+Mn(k)))'*invers_sigm(sm+1:sm+Mn(k),sm+1:sm+Mn(k))*Z_at(sm+1:sm+Mn(k),1);
       k=k+1;
    end
     
end
i=1;
SuM_AM=zeros(2*Nb,1);
while i<=Npmu
  
    SuM_AM=SuM_AM+SuM1(:,i);
     i=i+1;
end
V_hat_AM = (inv(G_AM))*SuM_AM;

%% step2 minimization with respect to gama -- estimation of theta

A= zeros(num_meas,2);

for i=1:1:Npmu
    
    sm=sum(Mn(1:i-1));
    
     o=Mn(i);
     j=1;
            while j<o
     
                 A(sm+j,1)= (H(sm+j,:))*V_hat_AM ;
                 A(sm+j,2)=-(H(sm+j+1,:))*V_hat_AM ;
                 
                 A(sm+j+1,1)=(H(sm+j+1,:))*V_hat_AM ;
                 A(sm+j+1,2)=(H(sm+j,:))*V_hat_AM ;
                j=j+2;
              
            end
    
end



    
    Gama_AM_est=zeros(2,1);
    Gama_takhmin=zeros(2,Npmu);
    k=1;
    for i=1:1:Nb
             if a(i)==1
         
             sm=sum(Mn(1:k-1));

                 
             T=(A(sm+1:sm+Mn(k),:))'*invers_sigm(sm+1:sm+Mn(k),sm+1:sm+Mn(k))*Z_at(sm+1:sm+Mn(k),1);
            j=1;
                     Gama_AM_est(:,1)=(1/(norm(T,2)))*(A(sm+1:sm+Mn(k),:))'*invers_sigm(sm+1:sm+Mn(k),sm+1:sm+Mn(k))*Z_at(sm+1:sm+Mn(k),1);
                     Gama_takhmin(1:2,k)=Gama_AM_est;
            j=j+1;
                    k=k+1;
               
             end
   
    end
    
 % gama
    gam11=Gama_takhmin(1,1);
    gam12=Gama_takhmin(2,1);
    gam21=Gama_takhmin(1,2);
    gam22=Gama_takhmin(2,2);
    gam31=Gama_takhmin(1,3);
    gam32=Gama_takhmin(2,3);
    gam41=Gama_takhmin(1,4);
    gam42=Gama_takhmin(2,4);
    gam51=Gama_takhmin(1,5);
    gam52=Gama_takhmin(2,5);
    gam61=Gama_takhmin(1,6);
    gam62=Gama_takhmin(2,6);
      
gama_new1=[gam11 -gam12;gam12 gam11];
gama_new2=[gam21 -gam22;gam22 gam21];
gama_new3=[gam31 -gam32;gam32 gam31];
gama_new4=[gam41 -gam42;gam42 gam41];
gama_new5=[gam51 -gam52;gam52 gam51];
gama_new6=[gam61 -gam62;gam62 gam61];



%% main approach for theta estimation --- approach with  EVD

   %EVD= ((A())'*invers_sigm()*A())
    
    
   % U()=(Q())'*(A())'*invers_sigm()*Z_at()
   
 %% objective function

  count=count+1
%    objective=zeros(Npmu,1)
for k=1:Npmu
     sm=sum(Mn(1:k-1));
       obj(d,k)= (Z_at(sm+1:sm+Mn(k),1)-A(sm+1:sm+Mn(k),:)*Gama_takhmin(1:2,k))'*invers_sigm(sm+1:sm+Mn(k),sm+1:sm+Mn(k))*(Z_at(sm+1:sm+Mn(k),1)-A(sm+1:sm+Mn(k),:)*Gama_takhmin(1:2,k));
%     if  objective(k) <=Tolerance
%          disp('its the end of the algorithm')
%          break
%    end
end
%   q=1;
i=1;
while i <= Npmu
    
    OBJ(d)=OBJ(d)+obj(d,i);
    i=i+1;
    
end
OBJ(d);
objective=abs(OBJ(d)-OBJ(d-1))/abs(OBJ(d));
if  objective <=Tolerance
    disp('------------------------------------------------------------------------------------------------')
         disp('its the end of the algorithm')
         break
end

    
  %% NEW Gama
 Gama_AM=zeros(num_meas,num_meas);
 k=1;
for i=1:1:Nb
     if a(i)==1
         
      sm=sum(Mn(1:k-1));

     o=Mn(k);
     j=1;
             if i==2
                        while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new1);
                       j=j+2;
                      
                        end
                        k=k+1;
              end
              
              
              
              
              if i==4
                         while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new2);
                          j=j+2;
                         
                         end
                         k=k+1;
              end
              
              if i==6
                         while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new3);
                          j=j+2;
                         
                         end
                         k=k+1;
              end
              
              if i==7
                         while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new4);
                          j=j+2;
                         
                         end
                         k=k+1;
              end
          
              if i==10
                         while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new5);
                          j=j+2;
                         
                         end
                         k=k+1;
              end
          
              
              if i==14
                         while j<=o
              
                          Gama_AM(sm+j:sm+j+1,sm+j:sm+j+1)= blkdiag(gama_new6);
                          j=j+2;
                         
                         end
                         k=k+1;
              end
          
              
                         
     end
    
end 

end


%%   final step
G1=[zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb),zeros(2*Nb,2*Nb)];

k=1;
for i=1:1:Nb
    if a(i)==1
        sm=sum(Mn(1:k-1));
        G1(:,(k-1)*2*Nb+1:(k-1)*2*Nb+2*Nb)=a(i)*(H(sm+1:sm+Mn(k),:))'*(Gama_AM(sm+1:sm+Mn(k),sm+1:sm+Mn(k)))'* invers_sigm(sm+1:sm+Mn(k),sm+1:sm+Mn(k))*(Gama_AM(sm+1:sm+Mn(k),sm+1:sm+Mn(k)))*H(sm+1:sm+Mn(k),:);
        k=k+1;
    end
end
i=1;
G_AM=zeros(2*Nb,2*Nb);
while i<=Npmu
  
    G_AM=G_AM+G1(:,(i-1)*2*Nb+1:(i-1)*2*Nb+2*Nb);
     i=i+1;
end

k=1;
 SuM1=zeros(2*Nb,2);
for i=1:1:Nb
    if a(i)==1
      sm=sum(Mn(1:k-1));
      SuM1(:,k)=a(i)*(H(sm+1:sm+Mn(k),:))'*(Gama_AM(sm+1:sm+Mn(k),sm+1:sm+Mn(k)))'*invers_sigm(sm+1:sm+Mn(k),sm+1:sm+Mn(k))*Z_at(sm+1:sm+Mn(k),1);
       k=k+1;
    end
     
end
i=1;
SuM_AM=zeros(2*Nb,1);
while i<=Npmu
  
    SuM_AM=SuM_AM+SuM1(:,i);
     i=i+1;
end
V_hat_AM = (inv(G_AM))*SuM_AM;



%%  result


disp('Gama_takhmin')
disp(Gama_takhmin)
disp('V_hat_AM ')
disp(V_hat_AM )
theta1=acos(Gama_takhmin(1,1))
theta2=acos(Gama_takhmin(1,2))
theta3=acos(Gama_takhmin(1,3))
theta4=acos(Gama_takhmin(1,4))
theta5=acos(Gama_takhmin(1,5))
theta6=acos(Gama_takhmin(1,6))
%% plot 
q1=figure('DefaultAxesFontName', font);
 k=1;
 plot(k,(180/pi)*delta_theta_spf1,'o','LineWidth',5,'MarkerSize',2,'color','g','LineWidth',5,'MarkerSize',2)
%  axis([0 14 0 50])
hold on
 plot(1,(180/pi)*theta1,'rs','LineWidth',5,'MarkerSize',2,'color','r','LineWidth',5,'MarkerSize',2)
%  xlabel('K')
 ylabel('Spoofing Angle(degree)')
 hold on

k=2;
 plot(k,(180/pi)*delta_theta_spf2,'o','LineWidth',5,'MarkerSize',2,'color','g','LineWidth',5,'MarkerSize',2)
%   axis([0 14 0 50])
hold on
 plot(2,(180/pi)*theta2,'rs','LineWidth',5,'MarkerSize',2,'color','r','LineWidth',5,'MarkerSize',2)
%  xlabel('K')
 ylabel('Spoofing Angle(degree)')
 hold on
 
k=3;
 plot(k,(180/pi)*delta_theta_spf3,'o','LineWidth',5,'MarkerSize',2,'color','g','LineWidth',5,'MarkerSize',2)
 hold on
 plot(3,(180/pi)*theta3,'rs','LineWidth',5,'MarkerSize',2,'color','r','LineWidth',5,'MarkerSize',2)
 xlabel('K')
 ylabel('Spoofing Angle(degree)')
 hold on


k=4;
 plot(k,(180/pi)*delta_theta_spf4,'o','LineWidth',5,'MarkerSize',2,'color','g','LineWidth',5,'MarkerSize',2)
 hold on
 plot(4,(180/pi)*theta4,'rs','LineWidth',5,'MarkerSize',2,'color','r','LineWidth',5,'MarkerSize',2)
%  xlabel('K')
 ylabel('Spoofing Angle(degree)')
 hold on
 

k=5;
 plot(k,(180/pi)*delta_theta_spf5,'o','LineWidth',5,'MarkerSize',2,'color','g','LineWidth',5,'MarkerSize',2)
hold on

plot(5,(180/pi)*theta5,'rs','LineWidth',5,'MarkerSize',2,'color','r','LineWidth',5,'MarkerSize',2)
%  xlabel('K')
 ylabel('Spoofing Angle(degree)')
 hold on
 
k=6;
 plot(k,(180/pi)*delta_theta_spf6,'o','LineWidth',5,'MarkerSize',2,'color','g','LineWidth',5,'MarkerSize',2)
 hold on
plot(6,(180/pi)*theta6,'rs','LineWidth',5,'MarkerSize',2,'color','r','LineWidth',5,'MarkerSize',2)
 xlabel('PMU number')
 ylabel('Spoofing Angle(degree)')
 hold on
  t=[1,2,3,4,5,6];
 set(gca,'XTick',t); 
 set(gca,'XTickLabel',t);
 axis([0 7 0 95])
hold on
legend('True value','Estimated value')
grid
print(q1,'-dpng','Figure1')



theta=[delta_theta_spf1;delta_theta_spf2;delta_theta_spf3;delta_theta_spf4;delta_theta_spf5;delta_theta_spf6];
teta_spf_estimation=[theta1;theta2;theta3;theta4;theta5;theta6];
estimation_error=  abs(teta_spf_estimation-theta)*(180/pi);

 disp('estimation error of spoofing theta (degree)=')
 disp(estimation_error)
est_error_percentage=zeros(6,1); 
 for i=1:6
 est_error_percentage(i)=abs((teta_spf_estimation(i)-theta(i))*(180/pi))/abs(theta(i)*(180/pi));
 end

 figure
plot(est_error_percentage)
 hold on
 
 
 MSE_1=zeros(6,1); 
 for i=1:6
 MSE_1(i)=((teta_spf_estimation(i)-theta(i))*(180/pi))^2;
 end
%  MSE=(1/6)*(sum(MSE_1))
 %%
error_AM=V-V_hat_AM;

X1=1:1:2*Nb ; 
Y1=V;
Y2=V_hat_AM;
X_length=length(X1);
figure
hold on
for nn=1:X_length
    x=X1(nn);
    y1=Y1(nn);
    y2=Y2(nn);
    plot(x,y1,'o','LineWidth',5,'MarkerSize',2,'color','g');
    hold on
    plot(x,y2,'pentagram','LineWidth',5,'MarkerSize',1,'color','r');
end
axis([0 2*Nb -2 2]);
xlabel('number of states')
ylabel('real part and imaginary part of each states')
title('state estimation (with attack detection)')
legend('true states','estimated states')


 hold on
 
q2=figure('DefaultAxesFontName', font);
axes;
hold on
x=1:1:2*Nb;
bar(x,[V,V_hat_AM])
xlabel('n')
% ylabel('V & V_hat_AM','fontsize',12)
xlabel('State number','fontsize',12)
ylabel('Real and imaginary parts of states','fontsize',12)
title('State estimation (with attack detection)','fontsize',12)
legend('True value','Estimated value','fontsize',12)
print(q2,'-dpng','Figure2')
teta_spf_estimation=[theta1;theta2;theta3;theta4;theta5;theta6]*(180/pi)
toc
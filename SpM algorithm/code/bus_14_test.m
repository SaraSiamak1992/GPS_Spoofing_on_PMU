clc
clear all
close all

  font = 'Times New Roman';
n=14;
Npmu=4;
%a=[0;1;0;0;0;1;1;0;1;0;0;0;0;0]; 

%% lines connected to bus n (bus with PMU)
n1=4 ;n2=4  ; n3=3;   n4=4  ; 

%% number of measurements
Mn=[(n1+1),(n2+1),(n3+1),(n4+1)];
num_meas=sum(Mn);

%% covariance matrix (sigma)

sg1=0.02 ;  sg2=0.03 ;sg3=0.03 ;sg4=0.02 ;  
Sg=[sg1,sg2,sg3,sg4];
sig=zeros(num_meas,1);

for k=1:1:Npmu
    sm=sum(Mn(1:k-1))
    sig(sm+1:sm+Mn(k))=Sg(k)*ones(Mn(k),1) ;
end

Sigma=diag(sig);
Ce=Sigma;
%% measurement noise
mu=0;
W=zeros(num_meas,1);
for i=1:1:Npmu
    sm=sum(Mn(1:i-1));
    W(sm+1:sm+Mn(i))=normrnd(mu,Sg(i),Mn(i),1);
end
%% measurements
H4pmu
S_states_14bus
% states
S;
S_real=real(S);
S_im=imag(S);
s_true=[S_real;S_im];
%S_correct=[real(S(1));imag(S(1));real(S(2));imag(S(2));real(S(3));imag(S(3));real(S(4));imag(S(4));real(S(5));imag(S(5));real(S(6));imag(S(6));real(S(7));imag(S(7))
%         real(S(8));imag(S(8));real(S(9));imag(S(9));real(S(10));imag(S(10));real(S(11));imag(S(11));real(S(12));imag(S(12));real(S(13));imag(S(13));real(S(14));imag(S(14))];

H;

m_true=H*S;
m=m_true;

%% Gama (true)  attack is occured on pmu p
p=input('enter number of  spoofed PMU (between 1 to 4)=') % spoofed PMU
Betta=input('spoofing theta (rad) on PMU  ; 0<theta<2pi = ')  
alfa=exp(Betta*j);

G_true=eye(num_meas,num_meas);
G=ones(num_meas,1);
 
for i=1:1:Npmu
    if i == p
    sm=sum(Mn(1:i-1));
    
    o=Mn(i);
    fi=1;
         while fi<=o
      
         G(sm+fi,1)=alfa ;
      
         fi=fi+1;
         end
   
    G_true=diag(G);
    end
end
G_true;

m_at=G_true*m_true;%+W;



%% state estimation before attack detection

% M=AS+e   A=Ar+jAm    S=E+jF           M=Mr+jMm=(Ar+jAm)(E+jF)+e
% M_tilda=[Mr;Mm]=[Ar -Am;Am  Ar][E;F]+e   A_tilda=[Ar -Am;Am  Ar]
% S_hat=[(A_tilda)'R(A_tilda)]^(-1)(A_tilda)'R[Mr;Mm]

 Hr=real(H);
 Him=imag(H);
 H_tilda=[Hr -Him
          Him  Hr];
      
    rem=real(m_at);
    imm=imag(m_at);
    Lm=[rem;imm];
    I= eye(2*num_meas);
    Cem=0.05*I;
    s_m=(pinv((ctranspose(H_tilda)*pinv(Cem)*H_tilda))*ctranspose(H_tilda)*pinv(Cem))*Lm;
    


 SS=[S_real
     S_im];
X1=1:1:2*n ; 
Y1=SS;
Y2=s_m;
X_length=length(X1);
w1=figure('DefaultAxesFontName', font);
axes;
hold on
for nn=1:X_length
    x=X1(nn);
    y1=Y1(nn);
    y2=Y2(nn);
    plot(x,y1,'o','LineWidth',5,'MarkerSize',2,'color','g');
    hold on
    plot(x,y2,'pentagram','LineWidth',5,'MarkerSize',1,'color','r');
end
axis([0 28 -2 2]);
xlabel('State number','fontsize',12)
ylabel('Real and imaginary parts of states','fontsize',12)
%title('State estimation (without attack detection)')
legend('True value','Estimated value','fontsize',12)
 hold on
 print(w1,'-dpng','Figure1')
 
 w1_1=figure('DefaultAxesFontName', font);
axes;
hold on
x=1:1:2*n;
bar(x,[SS,s_m])
xlabel('n','fontsize',12)
ylabel('Real and imaginary parts of states','fontsize',12)
xlabel('State number','fontsize',12)

%title('state estimation (without attack detection)')
legend('True value','Estimated value','fontsize',12)
grid on
print(w1_1,'-dpng','Figure1_1')
  
 Mr=real(m_at);
 Mm=imag(m_at);
 M_tilda=[Mr
          Mm];
 
tic
 %% initialization
ro=[0.61803;0.61803;0.61803;0.61803];

eps=10^(-5);

Thr=0.01;
I= eye(num_meas);
J=[zeros(100,Npmu),zeros(100,Npmu)];
teta_spf=zeros(Npmu,1);
a=zeros(100,Npmu);
b=zeros(100,Npmu);
%% Golden section search algorithm
for k=1:1:Npmu
   
X= ones(num_meas,1);
range_L=0;
range_R=2*pi;
  a(1,k)=range_L + (1-ro(k))*(range_R-range_L);
  b(1,k)=range_L + ro(k)*(range_R-range_L);
  
c=1;
      hghgf=0;
    while hghgf==0
%    while teta_spf(k)==0
   %% cumputation of J ----- { J(k,an) & J(k,bn) }
       
      %  J1  
       teta_spf_prime=a(c,k);
       
                       
                            sm=sum(Mn(1:k-1));
                           for q=1:1:Mn(k)
                           X(sm+q,1)= exp(1j*teta_spf_prime);
                           end
                     
    X;
    G=diag(X);
    G_prim=ctranspose(G);
    Y= I- H*((ctranspose(H)*Ce^(-1)*H)^(-1))*ctranspose(H)*Ce^(-1);
  
    r=Y*ctranspose(G)*m_at;
    J(c,k*2-1)=norm(r,2);
    J1=J(c,k*2-1);
  
        
      % J2
     X= ones(num_meas,1);
       teta_spf_prim=b(c,k);
      
             
                             sm=sum(Mn(1:k-1));
                           for q=1:1:Mn(k)
                           X(sm+q,1)= exp(1j*teta_spf_prim);
                           end  
             
    X;
    G=diag(X);
    Y= I- H*((ctranspose(H)*Ce^(-1)*H)^(-1))*ctranspose(H)*Ce^(-1);
    
    r=Y*ctranspose(G)*m_at;
    J(c,k*2)=norm(r,2);
    J2=J(c,k*2);
         
         
   %% updates
            if abs(range_L-range_R)< eps 
                 
   %% finding the spoofed theta of the Kth PMU
             X= ones(num_meas,1);
             J_prim_spf=zeros(2,1);
             D=[range_L;range_R];
             
                                   for z=1:1:2
                                      teta_s=D(z);

                                                            
                                                                       sm=sum(Mn(1:k-1));
                                                                       for q=1:1:Mn(k)
                                                                        X(sm+q,1)= exp(1j*teta_s);
                                                                       end
                                      G=diag(X) ;
                                      Y= I- H*((ctranspose(H)*Ce^(-1)*H)^(-1))*ctranspose(H)*Ce^(-1);
                                     
                                      r=Y*ctranspose(G)*m_at;
                                      J_prim_spf(z)=norm(r,2);
                                   end
             
            J_str=min(J_prim_spf(1),J_prim_spf(2)) ;
            o=find(J_prim_spf==J_str);
            teta_spf(k)=D(o);
             hghgf=1;
         
            else
                c=c+1
                     if  J1<J2
                         range_R = b(c-1,k);
                         b(c,k)  = a(c-1,k);
                         % J(k,bn+1)
                         J(c,k*2)= J(c-1,k*2-1);
                         
                         
                         a(c,k)=range_L + (1-ro(k))*(range_R-range_L);
                        
                     else
                         range_L   = a(c-1,k);
                         a(c,k)    = b(c-1,k);
                         % J(k,an+1)
                         J(c,k*2-1)= J(c-1,k*2);
                        
                         
                         b(c,k)=range_L + ro(k)*(range_R-range_L);
                        
                  
                     end
            end
        
        
        
    end
end

%% spofing theta 

                                 
 J_prim=zeros(Npmu,1);

 
 I= eye(num_meas);
 for k=1:1:Npmu
       X= ones(num_meas,1);
                          sm=sum(Mn(1:k-1));
                           for q=1:1:Mn(k)
                           X(sm+q,1)= exp(1j*teta_spf(k));
                           end
                           G=diag(X);
                            Y= I- H*((ctranspose(H)*Ce^(-1)*H)^(-1))*ctranspose(H)*Ce^(-1);
             
              r=Y*ctranspose(G)*m_at;
              J_prim(k)=norm(r,2);
 end 

            
    Jmin=min(J_prim)  ;
    o=find(J_prim==Jmin);
    teta_spf_estimation=teta_spf(o) ;                           
                                   
 X= ones(num_meas,1);
if teta_spf_estimation >= Thr
    disp('----------------------------------------------------------------')
    disp('Attack is occured!')
    teta_spf_guess=teta_spf
    disp('teta_spf_estimation (degree) =')
    disp(teta_spf_estimation*(180/pi))
    I= eye(num_meas);
  
     disp('a spoofing attack has occured on PMU :')
     disp(o)
     
 estimation_error=  abs(teta_spf_estimation-Betta)*(180/pi);

 disp('estimation error (degree)=')
 disp(estimation_error)

% plot
 w2=figure('DefaultAxesFontName', font);
 axes;
 plot(p,(180/pi)*Betta,'o','LineWidth',7,'MarkerSize',4,'color','r');%,'LineWidth',7,'MarkerSize',2)
 
 xlabel('PMU number','fontsize',12)
 ylabel(' \theta_s_p_f(degree)','fontsize',12)
 hold on
 plot(o,(180/pi)*teta_spf_estimation,'rs','LineWidth',7,'MarkerSize',3,'color','g');%,'LineWidth',7,'MarkerSize',4)
 
 t=[1,2,3,4];
 set(gca,'XTick',t); 
 set(gca,'XTickLabel',t);
  axis([0 5 , 0 50])
legend('\theta_s_p_f','Estimated \theta_s_p_f ','fontsize',12)
%title(' \theta_s_p_f estimation')
grid
  print(w2,'-dpng','Figure2')

                          sm=sum(Mn(1:o-1));
                           for q=1:1:Mn(o)
                           X(sm+q,1)= exp(1j*teta_spf_estimation);
                           end
                          
               
  
    X;
    G=diag(X);
    I= eye(2*num_meas,2*num_meas);
    m=ctranspose(G)*m_at;
    re=real(m);
    im=imag(m);
    L=[re;im];
    Ce=0.05*I;
    s=(pinv((ctranspose(H_tilda)*pinv(Ce)*H_tilda))*ctranspose(H_tilda)*pinv(Ce))*L; % pinv: pseudo inverse
    % s_true=((ctranspose(H)*Ce^(-1)*H)^(-1))*ctranspose(H)*Ce^(-1)*m_true
    
end
    
if teta_spf_estimation < Thr
   
disp('----------------------------------------------------------------')
    disp('condition is Normal')
     disp('teta_spf_estimation=')
    disp('zero')
    G=eye(num_meas);
  I= eye(2*num_meas,2*num_meas);
   m=ctranspose(G)*m_at;
    re=real(m);
    im=imag(m);
    L=[re;im];
    Ce=0.05*I;
    
    s=(pinv((ctranspose(H_tilda)*pinv(Ce)*H_tilda))*ctranspose(H_tilda)*pinv(Ce))*L;
     %s_true=((ctranspose(H)*Ce^(-1)*H)^(-1))*ctranspose(H)*Ce^(-1)*m_true
end

%      disp('s=')
%      disp(s')
     
     
     %% state estimation result
 Mr=real(m);
 Mm=imag(m);
 M_tilda=[Mr;Mm];
% S_correct;
%S_real_takhmin=(pinv((ctranspose(H_tilda)*pinv(Ce)*H_tilda))*ctranspose(H_tilda)*pinv(Ce))*M_tilda; 

%  figure 
%  x=1:1:28;
%   SS=[S_real
%      S_im];
%  plot(x,SS, x,s)
%      ylabel('state estimation (with attack detection)')
  SS=[S_real
     S_im];
X1=1:1:2*n ; 
Y1=SS;
Y2=s;
X_length=length(X1);
w3=figure('DefaultAxesFontName', font);
axes;
hold on
for nn=1:X_length
    x=X1(nn);
    y1=Y1(nn);
    y2=Y2(nn);
    plot(x,y1,'bd-','LineWidth',5,'MarkerSize',2,'color','g');
hold on
     plot(x,y2,'bd-','LineWidth',5,'MarkerSize',1,'color','b');

    xlabel('Number of states','fontsize',12)
    ylabel('Real and imaginary parts of states','fontsize',12)
    %title('State estimation (with attack detection)')
end
legend('True states','Estimated states','fontsize',12)
axis([0 2*n -2 2]);
 print(w3,'-dpng','Figure3')


w4=figure('DefaultAxesFontName', font);
axes;
hold on
x=1:1:2*n;
bar(x,[SS,s])
xlabel('n','fontsize',12)
ylabel('Real and imaginary parts of states','fontsize',12)
xlabel('State number','fontsize',12)

%title('state estimation (with attack detection)')
legend('True value','Estimated value','fontsize',12)
grid on
print(w4,'-dpng','Figure4')
toc 
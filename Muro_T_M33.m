clear
close all
clc

%% Diagrama de interacción - Muro T

%Input - Propiedades del elemento

fpc=210; %kg/cm2
betha_c=0.85; %Bloque de compresiones
fy=4200; %kg/cm2
Es=2.0000e+06; %kg/cm2
ey=0.0020; %Deformación unitaria del acero
ecu=0.003; %Según el ACI-318

%Input - Propiedades geométricas (cm)

b_1=675; %Ancho del ala
b_2=15; %Ancho del alma
h_1=15; %Altura del ala
h_2=415; %Altura total - Altura del alma h_2-h_1
Ag=b_1*h_1+b_2*(h_2-h_1); %cm2

%% Input - Distribución de aceros

n=36; %Número de capas de acero

%De arriba hacia abajo

As=[2.26; 2.26; 2.26; 2.26; 2.26;ones(12,1); 14.15; 14.15;...
     ones(12,1);2.26; 2.26; 2.26; 2.26; 2.26]; %Área de aceros - cm2

eid=[5; repmat(20,16,1); 10; 5; 10; repmat(20,16,1)]; %Distancia entre barras - cm

ei=[]; %cm

for i=1:length(eid)-1
    
    ei(1)=eid(1);
    ei=[ei ei(i)+eid(i+1)];
    
end

ei=ei';

n_c=1; %Avance del eje neutro

%% Compresión

%Cálculo del centroide de la sección bruta
%De arriba hacia abajo
y_c_1=b_1/2;
y_c_2=b_1-y_c_1;

c=0.00001:n_c:b_1+0.3*b_1; %Se puede cambiar según los puntos que se desea observar

e_i_a=[];

As_2=[];
ei_2=[];

As_1=As;

for i_3=0:length(As)-1

    As_2=[As_2; As(n-i_3)];
    ei_2=[ei_2;b_1-ei(n-i_3)];

end

flag=0;

%Carga máxima a compresión - Con estribos (ton)

P_n_max=0.8*(0.85*fpc*(Ag-sum(As))+fy*sum(As))/1000;
lim_1=(b_1-b_2)/2;
lim_2=(b_1+b_2)/2;

while flag==0 || flag==1
    
    a=zeros(length(c),1);
    P_n=zeros(length(c),1);
    P_real=zeros(length(c),1);
    M_n=zeros(length(c),1);
    M_real=zeros(length(c),1);
    phi_Mn=zeros(length(c),1);
    phi_Pn=zeros(length(c),1);

    for i=1:length(c)
        
        e_i=zeros(length(As),1);
        f_s=zeros(length(As),1);
        ms_i=zeros(length(As),1);

        a(i)=c(i)*betha_c;

        for i_1=1:length(As)

            %Cálculo de deformaciones en las barras

            e_i(i_1)=ecu*(ei(i_1)-c(i))/c(i);
            %Negativo - Compresión
            %Positivo - Tracción

            %Cálculo de fuerzas por capa de barras (kg)

            if abs(e_i(i_1))>ey

                if e_i(i_1)>0

                    f_s(i_1)=fy*As(i_1);

                else

                    f_s(i_1)=-fy*As(i_1);

                end

            else

                f_s(i_1)=Es*e_i(i_1)*As(i_1);

            end

            %Cálculo de fuerza de compresión del concreto (kg)

            if a(i)<lim_1

                P_con=0.85*fpc*a(i)*h_1;

            elseif lim_1<=a(i) && lim_2>a(i)

                P_con_1=0.85*fpc*h_1*(b_1-b_2)/2;
                P_con_2=0.85*fpc*(a(i)-(b_1-b_2)/2)*h_2;
                P_con=P_con_1+P_con_2;
                
            elseif a(i)>=lim_2
                
                P_con_1=0.85*fpc*b_2*(h_2-h_1);
                P_con_2=0.85*fpc*h_1*a(i);
                P_con=P_con_1+P_con_2;

            end

            %Momento generado por todas las capas de acero (kg-cm)

            ms_i(i_1)=f_s(i_1)*(y_c_1-ei(i_1));

        end

        %Almacenamiento de las deformaciones por cada eje neutro

        e_i_a=[e_i_a;e_i'];
        
        %Cálculo del momento nominal- Mn (ton.m)
        
        if a(i)<lim_1

            M_real(i)=(P_con*(y_c_1-a(i)/2)-sum(ms_i))/100000;

        elseif lim_1<=a(i) && lim_2>a(i)

            M_real(i)=(P_con_1*(y_c_1-(b_1-b_2)/4)+P_con_2*(y_c_1-(a(i)+(b_1-b_2)/2)/2)...
                -sum(ms_i))/100000;
            
        elseif a(i)>=lim_2
            
            M_real(i)=(P_con_2*(y_c_1-a(i)/2)-sum(ms_i))/100000;

        end

        %Fuerza axial total - Pn (ton)

        P_real(i)=(P_con-sum(f_s))/1000;
        %Negativo (Tracción)
        %Positivo (Compresión)
        
        if P_n_max>=P_real(i)
            
            P_n(i)=P_real(i);
            M_n(i)=M_real(i);
            
        else
            
            M_n(i)=0;
            P_n(i)=P_n_max;
            
        end

    end
    
    P_viga=0.1*fpc*Ag/1000; %kg
    M_b=max(M_n);
    p_1=find(M_n==M_b); %Posición del momento máximo
    P_b=P_n(p_1);
    
    P_tran_1=P_viga/0.7; %En el caso de muros con estribos
    P_tran_2=P_b/0.7;
    
    P_tran=min(P_tran_1,P_tran_2);
    
    %Cálculo de phi_Mn y phi_Pn
    
    for i_5=1:length(c)
        
        phi=0.9-0.2*P_n(i_5)/P_tran;
        
        if phi<=0.7
            
            phi=0.7;
            
        elseif phi>=0.9
            
            phi=0.9;
            
        end
        
        phi_Mn(i_5)=M_n(i_5)*phi;
        phi_Pn(i_5)=P_n(i_5)*phi;
        
    end
    
    %Almacenamiento de datos
    
    flag=flag+1;
    
    if flag==0 || flag==1
        
        M_n_neg=-M_n;
        P_n_neg=P_n;
        e_i_a_neg=e_i_a;
        
        phi_Mn_neg=-phi_Mn;
        phi_Pn_neg=phi_Pn;
        
        %Cambio para el eje negativo del diagrama
%         As=As_2;
%         ei=ei_2;
%         h_1_r=h_1;
%         h_2_r=h_2;
%         h_1=h_2_r;
%         h_2=h_1_r;
%         b_1_r=b_1;
%         b_2_r=b_2;
%         b_1=b_2_r;
%         b_2=b_1_r;
%         y_c_1=y_c_2;
        
    else
        
        M_n_pos=M_n;
        P_n_pos=P_n;
        e_i_a_pos=e_i_a;
        
        phi_Mn_pos=phi_Mn;
        phi_Pn_pos=phi_Pn;
        
    end
    

end

%Cambio de orden el gráfico para que se cierre el diagrama de interacción
M_n_pos_graf=[];
P_n_pos_graf=[];
phi_M_n_pos_graf=[];
phi_P_n_pos_graf=[];
n_3=length(M_n_pos);

for i_4=0:n_3-1
    
    M_n_pos_graf=[M_n_pos_graf; M_n_pos(n_3-i_4)];
    P_n_pos_graf=[P_n_pos_graf; P_n_pos(n_3-i_4)];
    phi_M_n_pos_graf=[phi_M_n_pos_graf; phi_Mn_pos(n_3-i_4)];
    phi_P_n_pos_graf=[phi_P_n_pos_graf; phi_Pn_pos(n_3-i_4)];
    
end

M_n_graf=[M_n_neg;M_n_pos_graf];
P_n_graf=[P_n_neg;P_n_pos_graf];

phi_M_n_graf=[phi_Mn_neg;phi_M_n_pos_graf];
phi_P_n_graf=[phi_Pn_neg;phi_P_n_pos_graf];

figure
plot(M_n_pos,P_n_pos,M_n_neg,P_n_neg)

figure
plot(M_n_graf,P_n_graf)
hold on
plot(phi_M_n_graf,phi_P_n_graf)

%% Solicitaciones del muro

load('T_M33.txt')
P_u_M33=T_M33(:,1); %ton
M_u_M33=T_M33(:,2); %ton.m

sz=15; %sz determina el tamano de los círculos

x_lim=2.5*10^3; %ton.m
y_lim_sup=3*10^3; %ton
y_lim_inf=1*10^3; %ton

figure
plot(M_n_graf,P_n_graf,'Color',[0 0.8 1],'LineWidth',2)
hold on
plot(phi_M_n_graf,phi_P_n_graf,'--','Color',[0 0.8 1],'LineWidth',2)
hold on
scatter(M_u_M33,P_u_M33,sz,'MarkerEdgeColor',[0.7 0.7 0.7],...
    'MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',0.01)
hold on
plot([-x_lim x_lim],[0 0],'k','LineWidth',0.9)
plot([0 0],[-y_lim_inf y_lim_sup],'k','LineWidth',0.9)
xticks(-x_lim:10^3:x_lim);
yticks(-y_lim_inf:1*10^3:y_lim_sup);
xlabel('Momento (ton-m)','FontSize',14,'FontName','Times New Roman');
ylabel('Carga Axial (ton)','FontSize',14,'FontName','Times New Roman','visible','on');
ytickformat('%.0f') %Se muestran cero decimales '0f' en el eje x
ax = gca;
ax.YAxis.Exponent = 3; %El exponente de los números del eje Y se ajustan a 3
ax.XAxis.Exponent = 3; %El exponente de los números del eje X se ajustan a 3
box off
set(gca,'linewidth',1.3,'FontSize',12,'FontName','Times New Roman')
axis([-x_lim x_lim -y_lim_inf y_lim_sup])

print('C_I_Conven_T_M33','-djpeg','-r300')

%Cálculo de la carga axial máxima

Pu_Nch=Ag*fpc*0.3/1000; %ton

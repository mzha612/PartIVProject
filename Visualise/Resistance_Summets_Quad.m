function Resistivity = Resistance_Summets_Quad(Vessel, cs)
%bc.bsD_p=[1 3 5 7 8 9]
% M=dlmread('input_charge_no_cell.txt'); % input parametre file
% id=str2num(job); % job number in the input file
% lambda=M(id,9);% Debey length coefficient
% lambda = asd;
% cs=M(id,12);%fix charge concentration
% cs = 5;
% if strcmp(flag,'shear')==1 || strcmp(flag,'shearsolid')==1
%     cs=x;
% end
% if strcmp(flag,'B')==1
%     cs=x;
% end
% ep=M(id,7);
% chi=M(id,16);
% chi = asd;
% perm=M(id,3);
% phi_f=M(id,6);
% V=M(id,10);
% gammap=M(id,13);
% gamman=M(id,14);
% ds=M(id,17);
% cstild=lambda*cs;
% lmbd=sqrt(2*lambda);
% h=1-ep;

%% Param

phi_f = 0.99; % Fluid phase fraction
mu_f = 10^-3; % viscosity of water, Pa.s
K = 10^10;
V = 1;
u0 = 1e-3;
H = Vessel.Radius;
D_Na = 1.62e-9;
D_Cl = 2.45e-9;
c0 = 927.388e23; % m^-3
esp = 678.2164e-12; % C(Vm)^-1;
T = 303.15; %K;
e = 1.6e-19; % C;
k = 1.38e-23; % JK^-1;

%% Non-Dim
h = Vessel.h;

perm = K*(H)^2/(phi_f * mu_f); % symbol chi
lambda = c0*e^2*H^2/(esp*k*T);
chi = (c0*k*T*H)/(mu_f*u0)/phi_f; % chi hat_l

% if lambda > 50000
%     disp('lambda large')
%     lambda = 50000;
% end
%
% if chi > 50000
%     disp('hartman large')
%     chi = 50000;
% end

gammap = (u0*H)/D_Na;
gamman = (u0*H)/D_Cl;

ds = 0;
cstild=lambda*cs;
lmbd=sqrt(2*lambda);

%%

%%
tol = 10^-10;
c_l=-cstild./(lmbd.^2).*sinh(lmbd-lmbd.*h)./sinh(lmbd);

c2_top=-c_l*sinh(lmbd.*h)./(2.*exp(lmbd).*sinh(lmbd-lmbd.*h));

c1_top=c2_top.*exp(2*lmbd);

if (isnan(c1_top) || isnan(c2_top) || isnan(c_l))
    c_l=-cs/2.*exp(-lmbd.*h);
    c2_top=-c_l/2*exp(2*lmbd*(h-1));
    c1_top=c2_top.*exp(2*lmbd);
    
    if abs(c_l) < tol
        c_l = 0;
    end
    
    if (isnan(c1_top))
        c1_top = 0;
    end    
end



if ((c1_top + c2_top + c_l) == 0) || isinf(exp(sqrt(2*lambda))) 
    Z1=0;
    Z2=0;
    Z3=-cs/2*(1-h);
    Z4=cs.^2/4*(1-h);
    Z5=0;
    Z7=+cs/(2*sqrt(perm))*(exp(-sqrt(perm))-exp(-h*sqrt(perm)));
    Z8=-cs/(2*sqrt(perm))*(exp(sqrt(perm))-exp(h*sqrt(perm)));
    
    
    W=(-h*perm*exp(sqrt(perm))+sqrt(perm)*exp(sqrt(perm)*h))./(h*perm*exp(-sqrt(perm))+sqrt(perm)*exp(-sqrt(perm)*h));
    tilQ=(-phi_f*chi/lambda-phi_f*2*chi/(perm-2*lambda))*pot_l(h)+phi_f*chi*c_l/lambda-phi_f*2*lambda*chi*cs/(perm*(perm-2*lambda));
    tilR=(-phi_f*chi/lambda-2*chi/(perm-2*lambda))*dpot_l(h);
    tilT=2*chi/(perm-2*lambda)*pot_l(1)+2*lambda*chi*cs/(perm*(perm-2*lambda));
    tilS=(tilR+tilT*h*perm)/(h*perm*exp(-sqrt(perm))+sqrt(perm)*exp(-sqrt(perm)*h));
    tilG=(h^2*perm/2*W*exp(-sqrt(perm))+h^2*perm/2*exp(sqrt(perm))-W*phi_f*(exp(-sqrt(perm)*h)-exp(-sqrt(perm)))-phi_f*(exp(sqrt(perm)*h)-exp(sqrt(perm))))^(-1);
    tilF=(-tilS*phi_f*(exp(-sqrt(perm)*h)-exp(-sqrt(perm)))-tilQ-tilT*perm*h.^2/2-phi_f*tilT+h.^2/2*perm*tilS*exp(-sqrt(perm)))*tilG;
    
    
else
    
    Z1=c_l./sqrt(2.*lambda).*sinh(h.*sqrt(2.*lambda));
    Z2=c_l.^2./2*(1/(2*sqrt(2*lambda))*sinh(h*2*sqrt(2*lambda))+h);
    Z3=1/sqrt(2*lambda)*(-c1_top*(exp(-sqrt(2*lambda))-exp(-sqrt(2*lambda)*h))+c2_top*(exp(sqrt(2*lambda))-exp(sqrt(2*lambda)*h)))-cs/2*(1-h);
    Z4=-c1_top.^2/(2*sqrt(2*lambda))*(exp(-2*sqrt(2*lambda))-exp(-2*sqrt(2*lambda)*h))+2*c1_top*c2_top*(1-h) ...
        +c2_top.^2/(2*sqrt(2*lambda))*(exp(2*sqrt(2*lambda))-exp(2*sqrt(2*lambda)*h)) ...
        +cs*c1_top/(sqrt(2*lambda))*(exp(-sqrt(2*lambda))-exp(-sqrt(2*lambda)*h)) ...
        -cs*c2_top/(sqrt(2*lambda))*(exp(sqrt(2*lambda))-exp(sqrt(2*lambda)*h))+cs.^2/4*(1-h);
    Z5=h.^2*Z1+1/lambda*Z1-h*c_l/lambda*cosh(h*sqrt(2*lambda));
    Z7=-c1_top/(sqrt(perm)+sqrt(2*lambda))*(exp(-sqrt(perm)-sqrt(2*lambda))-exp(-h*(sqrt(perm)+sqrt(2*lambda)))) ...
        +c2_top/(-sqrt(perm)+sqrt(2*lambda))*(exp(-sqrt(perm)+sqrt(2*lambda))-exp(h*(-sqrt(perm)+sqrt(2*lambda)))) ...
        +cs/(2*sqrt(perm))*(exp(-sqrt(perm))-exp(-h*sqrt(perm)));
    Z8=c1_top/(sqrt(perm)-sqrt(2*lambda))*(exp(sqrt(perm)-sqrt(2*lambda))-exp(h*(sqrt(perm)-sqrt(2*lambda)))) ...
        +c2_top/(sqrt(perm)+sqrt(2*lambda))*(exp(sqrt(perm)+sqrt(2*lambda))-exp(h*(sqrt(perm)+sqrt(2*lambda)))) ...
        -cs/(2*sqrt(perm))*(exp(sqrt(perm))-exp(h*sqrt(perm)));
    
    W=(-h*perm*exp(sqrt(perm))+sqrt(perm)*exp(sqrt(perm)*h))./(h*perm*exp(-sqrt(perm))+sqrt(perm)*exp(-sqrt(perm)*h));
    tilQ=(-phi_f*chi/lambda-phi_f*2*chi/(perm-2*lambda))*pot(h)+phi_f*chi*c_l/lambda-phi_f*2*lambda*chi*cs/(perm*(perm-2*lambda));
    tilR=(-phi_f*chi/lambda-2*chi/(perm-2*lambda))*dpot(h);
    tilT=2*chi/(perm-2*lambda)*pot(1)+2*lambda*chi*cs/(perm*(perm-2*lambda));
    tilS=(tilR+tilT*h*perm)/(h*perm*exp(-sqrt(perm))+sqrt(perm)*exp(-sqrt(perm)*h));
    tilG=(h^2*perm/2*W*exp(-sqrt(perm))+h^2*perm/2*exp(sqrt(perm))-W*phi_f*(exp(-sqrt(perm)*h)-exp(-sqrt(perm)))-phi_f*(exp(sqrt(perm)*h)-exp(sqrt(perm))))^(-1);
    tilF=(-tilS*phi_f*(exp(-sqrt(perm)*h)-exp(-sqrt(perm)))-tilQ-tilT*perm*h.^2/2-phi_f*tilT+h.^2/2*perm*tilS*exp(-sqrt(perm)))*tilG;
    
    
end

    function f=pot(y)
        f=c1_top*exp(-lmbd*y)+c2_top*exp(lmbd*y)-cstild/(lmbd.^2);
    end

    function f=pot_l(y)
        f=-cstild/(lmbd.^2);
    end

    function f=dpot(y)
        f=-lmbd*c1_top*exp(-lmbd*y)+lmbd*c2_top*exp(lmbd*y);
    end

    function f=dpot_l(y)
        f=0;
    end


chek=0;
if(gammap==0 || gamman==0)
    chek=1;
end
switch chek
    case 0
        if ds==0
            
            B=(-2*V*Z1+perm*(Z5-2*Z3/perm)*(exp(-sqrt(perm))*V*W*tilG+exp(sqrt(perm))*V*tilG)+2*Z7*V*W*tilG+2*Z8*V*tilG) ...
                *(perm*(Z5-2*Z3/perm)*(exp(-sqrt(perm))*(W*tilF-tilS)+exp(sqrt(perm))*tilF+tilT)-2*phi_f*chi/lambda*Z2 ...
                +phi_f*chi*c_l*Z1+2*Z7*(W*tilF-tilS)+2*Z8*tilF+4*chi/(perm-2*lambda)*Z4+Z3*4*lambda*chi*cs/(perm*(perm-2*lambda)) ...
                +(1/gammap+1/gamman)*(2+(Z2+Z4)/2)+(1/gamman-1/gammap)*(Z1+Z3) ).^(-1);
        else
            B=0;
        end
    case 1
        gammap=1.7;
        gamman=2.5
        B=(-2*V*Z1+perm*(Z5-2*Z3/perm)*(exp(-sqrt(perm))*V*W*tilG+exp(sqrt(perm))*V*tilG)+2*Z7*V*W*tilG+2*Z8*V*tilG) ...
            *(perm*(Z5-2*Z3/perm)*(exp(-sqrt(perm))*(W*tilF-tilS)+exp(sqrt(perm))*tilF+tilT)-2*phi_f*chi/lambda*Z2 ...
            +phi_f*chi*c_l*Z1+2*Z7*(W*tilF-tilS)+2*Z8*tilF+4*chi/(perm-2*lambda)*Z4+Z3*4*lambda*chi*cs/(perm*(perm-2*lambda)) ...
            +(1/gammap+1/gamman)*(2+(Z2+Z4)/2)+(1/gamman-1/gammap)*(Z1+Z3) ).^(-1);
        % B=1;
end
if ((c1_top + c2_top + c_l) == 0) || isinf(exp(sqrt(2*lambda)))
    Q=(-phi_f*chi*B/lambda-phi_f*2*chi*B/(perm-2*lambda))*pot_l(h)+V+phi_f*chi*B*c_l/lambda-phi_f*2*lambda*chi*B*cs/(perm*(perm-2*lambda));
    R=(-phi_f*chi*B/lambda-2*chi*B/(perm-2*lambda))*dpot_l(h);
    T=2*chi*B/(perm-2*lambda)*pot_l(1)+2*lambda*chi*B*cs/(perm*(perm-2*lambda));
    %W=(-h*perm*exp(sqrt(perm))+sqrt(perm)*exp(sqrt(perm)*h))./(h*perm*exp(-sqrt(perm))+sqrt(perm)*exp(-sqrt(perm)*h));
    S=(R+T*h*perm)/(h*perm*exp(-sqrt(perm))+sqrt(perm)*exp(-sqrt(perm)*h));
    D2=(-S*phi_f*(exp(-sqrt(perm)*h)-exp(-sqrt(perm)))-Q-T*perm*h.^2/2-phi_f*T+h.^2/2*perm*S*exp(-sqrt(perm)))*(h^2*perm/2*W*exp(-sqrt(perm))+h^2*perm/2*exp(sqrt(perm))-W*phi_f*(exp(-sqrt(perm)*h)-exp(-sqrt(perm)))-phi_f*(exp(sqrt(perm)*h)-exp(sqrt(perm))))^(-1);
    D1=W*D2-S;
    A=perm*(D1.*exp(-sqrt(perm))+D2.*exp(sqrt(perm))+T);
else
    Q=(-phi_f*chi*B/lambda-phi_f*2*chi*B/(perm-2*lambda))*pot(h)+V+phi_f*chi*B*c_l/lambda-phi_f*2*lambda*chi*B*cs/(perm*(perm-2*lambda));
    R=(-phi_f*chi*B/lambda-2*chi*B/(perm-2*lambda))*dpot(h);
    T=2*chi*B/(perm-2*lambda)*pot(1)+2*lambda*chi*B*cs/(perm*(perm-2*lambda));
    %W=(-h*perm*exp(sqrt(perm))+sqrt(perm)*exp(sqrt(perm)*h))./(h*perm*exp(-sqrt(perm))+sqrt(perm)*exp(-sqrt(perm)*h));
    S=(R+T*h*perm)/(h*perm*exp(-sqrt(perm))+sqrt(perm)*exp(-sqrt(perm)*h));
    D2=(-S*phi_f*(exp(-sqrt(perm)*h)-exp(-sqrt(perm)))-Q-T*perm*h.^2/2-phi_f*T+h.^2/2*perm*S*exp(-sqrt(perm)))*(h^2*perm/2*W*exp(-sqrt(perm))+h^2*perm/2*exp(sqrt(perm))-W*phi_f*(exp(-sqrt(perm)*h)-exp(-sqrt(perm)))-phi_f*(exp(sqrt(perm)*h)-exp(sqrt(perm))))^(-1);
    D1=W*D2-S;
    A=perm*(D1.*exp(-sqrt(perm))+D2.*exp(sqrt(perm))+T);
    
    
end

    function f=vlum(y)
        f=A.*y.^2./2+phi_f.*chi*B.*c_l./lambda.*(1-cosh(sqrt(2.*lambda).*y))+V;
    end
    function f=vEGL(y)
        f=D1*exp(-sqrt(perm).*y)+D2*exp(sqrt(perm).*y)+2*chi*B/(perm-2*lambda).*pot(y)-A/perm+2*lambda*chi*B*cs./(perm*(perm-2*lambda));
        
    end
    function f=vEGL_bot(y)
        
        %         f=D2*exp(-sqrt(perm).*y)+D1*exp(sqrt(perm).*y)+2*chi*B/(perm-2*lambda).*pot_bot(0.*x,y)-A/perm+2*lambda*chi*B*cs./(perm*(perm-2*lambda));
        f=D2*exp(-sqrt(perm).*y)+D1*exp(sqrt(perm).*y)+2*chi*B/(perm-2*lambda).*pot_bot(0,y)-A/perm+2*lambda*chi*B*cs./(perm*(perm-2*lambda));
        
    end
%     function f=pot_top(x,y)
%         f=c1_top*exp(-lmbd*y)+c2_top*exp(lmbd*y)-cstild/(lmbd.^2)+B*x;
%     end
    function f=pot_bot(x,y)
        f=c2_top*exp(-lmbd*y)+c1_top*exp(lmbd*y)-cstild/(lmbd.^2)+B*x;
    end
%     function f=pot_mid(x,y)
%         f=c_l*cosh(lmbd*y)+B*x;
%     end
%
%     function f=Ilum(y)
%         f=-2.*vlum(y).*pot_mid(0,y)-B.*(1./gammap.*(1-pot_mid(0,y)+pot_mid(0,y).*pot_mid(0,y)./2)+1./gamman.*(1+pot_mid(0,y)+pot_mid(0,y).*pot_mid(0,y)./2)).*(y-y+1);
%     end
%     function f=IEGL(y)
%         f=-2.*vEGL(y).*pot_top(0,y)-B.*(1./gammap.*(1-pot_top(0,y)+pot_top(0,y).*pot_top(0,y)./2)+1./gamman.*(1+pot_top(0,y)+pot_top(0,y).*pot_top(0,y)./2)).*(y-y+1);
%
%     end
%     function f=IEGL_bot(y)
%
%         f=-2.*vEGL_bot(y).*pot_bot(0,y)-B.*(1/gammap*(1-pot_bot(0,y)+pot_bot(0,y).*pot_bot(0,y)./2)+1/gamman.*(1+pot_bot(0,y)+pot_bot(0,y).*pot_bot(0,y)./2)).*(y-y+1);
%
%     end

% %___________solid phase solution_______________
%     function f=dvEGL(y)
%         f=-sqrt(perm).*D1*exp(-sqrt(perm).*y)+sqrt(perm).*D2*exp(sqrt(perm).*y)+2*chi*B/(perm-2*lambda).*dpot(y);
%
%     end
% M=((1-phi_f)./phi_f+1).*dvEGL(h)+chi*B/lambda*dpot(h)-((1-phi_f)./phi_f+1)*A*h;
% N=chi*B/lambda.*pot_top(0,1)-((1-phi_f)./phi_f+1)*A/2-M;
%     function f=u_top(y)
%         f=-vEGL(y)*((1-phi_f)./phi_f+1)-chi*B/lambda.*pot_top(0,y)+((1-phi_f)./phi_f+1).*A.*y.^2./2+M.*y+N;
%
%     end


%_____________________________________________
% switch flag
%     case 'vmid'
%         z=vlum(y);
%     case 'vtop'
%         z=vEGL(y);
%     case 'utop'
%         z=u_top(y);
%     case 'vbot'
%         z=vEGL_bot(y);
% %     case 'phi_top'
% %         z=pot_top(x,y);
% %     case 'cp_top'
% %         z=1-pot_top(x,y)+pot_top(x,y).^2/2;
% %     case 'cn_top'
% %         z=1+pot_top(x,y)+pot_top(x,y).^2/2;
% %     case 'phi_mid'
% %         z=pot_mid(x,y);
% %     case 'cp_mid'
% %         z=1-pot_mid(x,y)+pot_mid(x,y).^2/2;
% %     case 'cn_mid'
% %         z=1+pot_mid(x,y)+pot_mid(x,y).^2/2;
% %     case 'phi_bot'
% %         z=pot_bot(x,y);
% %     case 'cp_bot'
% %         z=1-pot_bot(x,y)+pot_bot(x,y).^2/2;
% %     case 'cn_bot'
% %         z=1+pot_bot(x,y)+pot_bot(x,y).^2/2;
%     case 'Imid'
%         z=Ilum(y);
%     case 'Itop'
%         z=IEGL(y);
%     case 'Ibot'
%         z=IEGL_bot(y);
%     case 'shear'
%         %shear_stress=-(-sqrt(perm).*D1*exp(-sqrt(perm))+sqrt(perm)*D2*exp(sqrt(perm))+2*chi*B/(perm-2*lambda).*(-lmbd*c1_top*exp(-lmbd)+lmbd*c2_top*exp(lmbd)))
%         z=-(-sqrt(perm).*D1*exp(-sqrt(perm))+sqrt(perm)*D2*exp(sqrt(perm))+2*chi*B/(perm-2*lambda).*(-lmbd*c1_top*exp(-lmbd)+lmbd*c2_top*exp(lmbd)));
%     case 'shearsolid'
%         z=-(-dvEGL(1)*((1-phi_f)./phi_f+1)-chi*B/lambda.*dpot(1)+((1-phi_f)./phi_f+1).*A+M);
%     case 'B'
%         z=B;
%
%
% end

    function f=vlum_l(y)
        f=A.*y.^2./2+V;
    end
    function f=vEGL_l(y)
        f=D1*exp(-sqrt(perm).*y)+D2*exp(sqrt(perm).*y)+2*chi*B/(perm-2*lambda).*pot_l(y)-A/perm+2*lambda*chi*B*cs./(perm*(perm-2*lambda));
        
    end
    function f=vEGL_bot_l(y)
        
        %         f=D2*exp(-sqrt(perm).*y)+D1*exp(sqrt(perm).*y)+2*chi*B/(perm-2*lambda).*pot_bot(0.*x,y)-A/perm+2*lambda*chi*B*cs./(perm*(perm-2*lambda));
        f=D2*exp(-sqrt(perm).*y)+D1*exp(sqrt(perm).*y)+2*chi*B/(perm-2*lambda).*pot_bot_l(0,y)-A/perm+2*lambda*chi*B*cs./(perm*(perm-2*lambda));
        
    end
    function f=pot_bot_l(x,y)
        f=0-cstild/(lmbd.^2)+B*x;
    end

if ((c1_top + c2_top + c_l) == 0) || isinf(exp(sqrt(2*lambda)))
    step = (1-h)/10;
    ytop = h:step:1;
    ymid = -h:step:h;
    ybot = -1:step:-h;
    
    v_top = vEGL_l(ytop);
    v_mid = vlum_l(ymid);
    v_bot = vEGL_bot_l(ybot);
    
else
    step = (1-h)/20;
    ytop = h:step:1;
    ymid = -h:step:h;
    ybot = -1:step:-h;
    disp('hi')
    v_top = vEGL(ytop);
    v_mid = vlum(ymid);
    v_bot = vEGL_bot(ybot);
end
%     function f=vEGL2(y)
%         f=D1*exp(-sqrt(perm).*y)+D2*exp(sqrt(perm).*y)+2*chi*B/(perm-2*lambda).*(c1_top*exp(-lmbd*y)+c2_top*exp(lmbd*y)-cstild/(lmbd.^2))-A/perm+2*lambda*chi*B*cs./(perm*(perm-2*lambda));
%
%     end
% function f=pot(y)
%      f=;
%     end
if isempty(v_top)
    F1 = 0;
else
    F1 = trapz(ytop,v_top);
end

if isempty(v_mid)
    F2 = 0;
else
    
    F2 = trapz(ymid,v_mid);
end

if isempty(v_bot)
    F3 = 0;
else
    F3 = trapz(ybot,v_bot);
%     F1 = trapz(ytop,v_top);
end


Flowrate = F1 + F2 + F3;
% Flowrate = (integral(vEGL2,h,1) + integral(vlum,-h,h) + integral(vEGL_bot,-1,-h));

%% Calculate resistance
% Division of the pressure gradient A, with the flow rate
Resistivity = -A/Flowrate;

end
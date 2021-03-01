function [PA_new] = meancurvature(i,ico,hemi,v_A,area,path,TR_A,F_A)

nsample=4^ico*10+2; %number of particles

[v_B,f_B]=read_vtk(strcat('/',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.white.s5.vtk'));
F_B=f_B+1;
FV.vertices = v_B;
FV.faces = F_B;
F_B=double(F_B);
TR_B=triangulation(F_B,v_B);
[Cmean_B, Cgaussian_B, Dir1_B, Dir2_B, Lambda1_B, Lambda2_B] = patchcurvature(FV,false);
Cmean_B = zscore(Cmean_B);
Cmean_B(Cmean_B < -3) = -3 - (1 - exp(3 + Cmean_B(Cmean_B < -3)));
Cmean_B(Cmean_B > 3) = 3 + (1 - exp(3 - Cmean_B(Cmean_B > 3)));
Cmean_B = Cmean_B / std(Cmean_B);
Cmean_B(Cmean_B < -3) = -3 - (1 - exp(3 + Cmean_B(Cmean_B < -3)));
Cmean_B(Cmean_B > 3) = 3 + (1 - exp(3 - Cmean_B(Cmean_B > 3)));
% Cmean_B = Cmean_B*10;

FV.vertices = v_A;
FV.faces = F_A;

[Cmean_A, Cgaussian_A, Dir1_A, Dir2_A, Lambda1_A, Lambda2_A] = patchcurvature(FV,false);
Cmean_A = zscore(Cmean_A);
Cmean_A(Cmean_A < -3) = -3 - (1 - exp(3 + Cmean_A(Cmean_A < -3)));
Cmean_A(Cmean_A > 3) = 3 + (1 - exp(3 - Cmean_A(Cmean_A > 3)));
Cmean_A = Cmean_A / std(Cmean_A);
Cmean_A(Cmean_A < -3) = -3 - (1 - exp(3 + Cmean_A(Cmean_A < -3)));
Cmean_A(Cmean_A > 3) = 3 + (1 - exp(3 - Cmean_A(Cmean_A > 3)));
% Cmean_A = Cmean_A*10;

[v_S,f_S]=read_vtk(strcat('/',path,'/oursphere/xyz/',num2str(i),'/',num2str(i),'ico',num2str(ico),'.real.vtk'));
F_S=f_S+1;
% compute the bary for after surgery (used in map_A, PA and PA_new)
cmd1 = strcat('/',path,'/SphericalRemesh -s /',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.sphere.vtk -r /',path,'/oursphere/xyz/',num2str(i),'/',num2str(i),'ico',num2str(ico),'.real.vtk --bary /',path,'/bary/before/',num2str(i),'/',num2str(i),'curv_baryico',num2str(ico),'real.txt');
cmd2 = strcat('/',path,'/SphericalRemesh -s /',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.sphere.vtk -r /',path,'/oursphere/ico',num2str(ico),'.vtk --bary /',path,'/bary/after/',num2str(i),'/',num2str(i),'curv_baryico',num2str(ico),'real.txt');
system(cmd1);
system(cmd2);

neighbor = cell(length(v_S), 1);
for j = 1: size(F_S, 1)
    neighbor{F_S(j,1)} = [neighbor{F_S(j,1)}, F_S(j, 2: 3)];
    neighbor{F_S(j,2)} = [neighbor{F_S(j,2)}, F_S(j, 1:2:3)];
    neighbor{F_S(j,3)} = [neighbor{F_S(j,3)}, F_S(j, 1: 2)];
end
for j = 1: length(v_S)
    neighbor{j} = unique(neighbor{j});
    nn(j)=length(neighbor{j});
end
PB=load(strcat('/',path,'/bary/after/',num2str(i),'/',num2str(i),'curv_baryico',num2str(ico),'real.txt'));
PB(:,1)=PB(:,1)+1;
PA=load(strcat('/',path,'/bary/before/',num2str(i),'/',num2str(i),'curv_baryico',num2str(ico),'real.txt'));
PA(:,1)=PA(:,1)+1;
map_B=zeros(size(PB,1),4);
map_A=zeros(size(PA,1),4);
for j=1:size(PB,1)    
    BC_B=PB(j,2:4);
    ID_B=PB(j,1);
    % update features
    map_B(j,:) = [BC_B*v_B(F_B(ID_B,:),:) BC_B*Cmean_B(F_B(ID_B,:),:)];
end
[~,d] = dsearchn(v_A,map_B(:,1:3));
R=1./(1+exp(-(7-d)));
    % R=ones(size(PA,1),1);

%% gradient descend
k=1;
eta = 1;
alpha=0.05;
PA_new=PA;
map_A_prev=map_A;
niter=10;
energy = zeros(niter,1);
energy_sim = zeros(niter,1);
energy_reg = zeros(niter,1);

while k <= niter %&& alpha > 1e-6
    % update features
    tic;
    for j=1:nsample
        BC_A=PA_new(j,2:4);
        ID_A=PA_new(j,1);
        % update features
        map_A(j,:) = [BC_A*v_A(F_A(ID_A,:),:) BC_A*Cmean_A(F_A(ID_A,:),:)];
    end
    ticFeature = toc;
    x_new = map_A;

    if k == 1
        x0 = map_A(:,1:3);
    end

    % energy function
    tic;
    energy(k) = 0;
    energy_sim(k) = 0;
    energy_reg(k) = 0;
    for j=1:nsample
        e_sim=norm(x_new(j,1:3)-map_B(j,1:3))^2;
        e_sim=e_sim*R(j);
        e_reg=0;
        coordn_B=map_B(neighbor{j},1:3);
        coordn_A=map_A(neighbor{j},1:3);

        tmp_reg=vecnorm(x_new(j,1:3)-coordn_A(:,1:3),2,2)-vecnorm(map_B(j,1:3)-coordn_B(:,1:3),2,2);
        tmp_reg=dot(tmp_reg,tmp_reg);
        e_reg=e_reg+tmp_reg;
        e_reg=e_reg/nn(j);

        energy_sim(k)=energy_sim(k)+e_sim;
        energy_reg(k)=energy_reg(k)+e_reg;
    end
    energy(k) = energy(k) + energy_sim(k) +eta*energy_reg(k);
    ticEnergy = toc;

    if k>1 && energy(k) > energy(k-1)
        % update features
        map_A=map_A_prev;
        alpha=alpha/10;
        x_new = map_A;
%         k=k-1;
    else
        map_A_prev=map_A;
        PA=PA_new;
    end
    %see energy
    fprintf('[%0#3d] energy for H: %.4f (%.4f + 0.5 * %.4f), alpha=%.4f',k,energy(k),energy_sim(k),energy_reg(k),alpha);

    % debug
    tic;
    for j=1:nsample
        tic;
        ID_A=PA(j,1);
        BC_A=PA(j,2:4);

        v_Acoord1=v_A(F_A(ID_A,1),:);        
        v_Acoord2=v_A(F_A(ID_A,2),:);
        v_Acoord3=v_A(F_A(ID_A,3),:);
        v_Acur1 = Cmean_A(F_A(ID_A,1),:);
        v_Acur2 = Cmean_A(F_A(ID_A,2),:);
        v_Acur3 = Cmean_A(F_A(ID_A,3),:);

        u=map_A(j,:);
        n_t=TR_A.faceNormal(ID_A);

        g1=cross(n_t,v_Acoord3-v_Acoord2)/2;
        g2=cross(n_t,v_Acoord1-v_Acoord3)/2;
        g3=cross(n_t,v_Acoord2-v_Acoord1)/2;

        G = [];
        G1 = [];
        for q = 1:3
            G = [v_Acoord1(q), v_Acoord2(q), v_Acoord3(q)]*[g1;g2;g3]/area(ID_A);
            G1 = [G1; G];
        end  
        p = 1;
        G2 = [v_Acur1(p), v_Acur2(p), v_Acur3(p)]*[g1;g2;g3]/area(ID_A);

        g_sim=(u(:,1:4)-map_B(j,1:4))*[G1;G2];
        g_sim=g_sim*2*R(j);


        coordn_B=map_B(neighbor{j},:);
        coordn_A=map_A(neighbor{j},:);

        term00=coordn_A(:,1:3)-u(:,1:3);
        term01=vecnorm(term00,2,2);
        term1=term01 - vecnorm(map_B(j,1:3)-coordn_B(:,1:3),2,2);
        term2=1./term01;
        term3=term00*G1;
        g_reg=(term1.*term2)'*term3/nn(j);

        g_total=-alpha*(g_sim+eta*g_reg);
        tic;
        [new_ID,new_BC]=PM(ID_A,BC_A,g_total,TR_A,v_A,F_A);
        PA_new(j,:)=[new_ID new_BC];
    end
     ticGrad = toc;

    fprintf('(%.2f = %.2f, %.2f, %.2f)\n', ticFeature+ticEnergy+ticGrad, ticFeature,ticEnergy,ticGrad);

    k=k+1;
end

    % plot(1:k-1,energy,'-*');

end


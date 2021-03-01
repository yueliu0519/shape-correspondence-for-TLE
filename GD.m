function [PA_new] = GD(niter, eta, alpha, nsample, i, ico, hemi,featype,area, map_B,R,path,TR_A,v_A,F_A)
%load initial PA of after surgery
PA=load(strcat('/',path,'/bary/before/',num2str(i),'/',num2str(i),'baryico',num2str(ico),'real.txt'));
PA(:,1)=PA(:,1)+1;
map_A=zeros(size(PA,1),3);
PA_new=PA;
map_A_prev=map_A;
% initialize energy
energy = zeros(niter,1);
energy_sim = zeros(niter,1);
energy_reg = zeros(niter,1);
% neighborhood info
if ico == 0
    [v_S,f_S]=read_vtk(strcat('/',path,'/oursphere/ico',num2str(ico),'.vtk'));
    F_S=f_S+1;
else
    [v_S,f_S]=read_vtk(strcat('/',path,'/oursphere/',featype,'/',num2str(i),'/',num2str(i),'ico',num2str(ico),'.real.vtk'));
    F_S=f_S+1;
end
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
k=1;

while k <= niter %&& alpha > 1e-6
    % update features
    tic;
    for j=1:nsample
        BC_A=PA_new(j,2:4);
        ID_A=PA_new(j,1);
        % update features
        map_A(j,:)=BC_A*v_A(F_A(ID_A,:),:);
    end
    ticFeature = toc;
    x_new = map_A;
    if k == 1
        x0 = x_new;
    end

% energy function
    tic;
    energy(k) = 0;
    energy_sim(k) = 0;
    energy_reg(k) = 0;
    for j=1:nsample
        e_sim=norm(x_new(j,:)-map_B(j,:))^2;
        e_sim=e_sim*R(j);
        e_reg=0;
        coordn_B=map_B(neighbor{j},:);
        coordn_A=map_A(neighbor{j},:);

        tmp_reg=vecnorm(x_new(j,:)-coordn_A,2,2)-vecnorm(map_B(j,:)-coordn_B,2,2);
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
        k=k-1;
    else
        map_A_prev=map_A;
        PA=PA_new;
    end
%see energy
    fprintf('[%0#3d] energy for xyz: %.4f (%.4f + 0.5 * %.4f), alpha=%.4f',k,energy(k),energy_sim(k),energy_reg(k),alpha);

% debug
    tic;
    parfor j=1:nsample
        tic;
        ID_A=PA(j,1);
        BC_A=PA(j,2:4);

        v_Acoord1=v_A(F_A(ID_A,1),:);        
        v_Acoord2=v_A(F_A(ID_A,2),:);
        v_Acoord3=v_A(F_A(ID_A,3),:);

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

        g_sim=(u-map_B(j,:))*G1;
        g_sim=g_sim*2*R(j);


        coordn_B=map_B(neighbor{j},:);
        coordn_A=map_A(neighbor{j},:);

        term00=coordn_A-u;
        term01=vecnorm(term00,2,2);
        term1=term01 - vecnorm(map_B(j,:)-coordn_B,2,2);
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

end
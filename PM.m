function [new_ID,new_BC]=PM(j,BC,g,TR,v,F)
if norm(g)<1e-15
    new_ID=j;
    new_BC=BC;
%     X=(BC*v(F(j,:),:));
    return
end
%% compute normal n_t of t
n_t=TR.faceNormal(j);
%% compute g x n for the plane normal n_p
n_p=cross(n_t,g);
n_p=n_p/norm(n_p);
%% define the plane P
x1=(BC*v(F(j,:),:));
d=-n_p * x1';
a=n_p(1);
b=n_p(2);
c=n_p(3);
%% intersection
v1=v(F(j,1),:);
v2=v(F(j,2),:);
v3=v(F(j,3),:);
[X,Y,x_b,y_b,vertex_x,vertex_y,edge_x,edge_y]=intersection(v1,v2,v3,a,b,c,d);
%% find positive sign
if (X-x1)*g' < 0
    X=Y;
    edge_x=edge_y;
    vertex_x=vertex_y;
    x_b=y_b;
end
%% ||g||
leng = norm(g);
len = norm(X-x1);

if len > leng
    X = (X - x1) * leng / len + x1;
    len = leng;
end
%% recursive function
new_ID = j;
exc_ID = new_ID;

while len < leng
    if x_b==1 % 4-1
        t_ID=edgeAttachments(TR,F(new_ID,edge_x(1)),F(new_ID,edge_x(2)));
        new_ID = setdiff(t_ID{1,1},new_ID);
%         new_ID = full(edgeA(F(new_ID,edge_x(2)),F(new_ID,edge_x(1))));
    else % 4-2
        X_ID=F(new_ID,vertex_x);
        t_ID=setdiff(cell2mat(vertexAttachments(TR,X_ID)),exc_ID);
        for k=1:length(t_ID)
            if t_ID(k) == new_ID
                continue;
            end
            XC=setdiff(F(t_ID(k),:),X_ID);
            s = [v(XC,:) ones(2,1)] * [a b c d]';
            if s(1) * s(2) <= 0
                new_ID = t_ID(k);
                break;
            end
        end
        exc_ID = t_ID;
    end
    v1=v(F(new_ID,1),:);
    v2=v(F(new_ID,2),:);
    v3=v(F(new_ID,3),:);
    x1 = X;
    [X,Y,x_b,y_b,vertex_x,vertex_y,edge_x,edge_y]=intersection(v1,v2,v3,a,b,c,d);
    if norm(x1-X) < norm(x1-Y)
        X=Y;
        edge_x=edge_y;
        vertex_x=vertex_y;
        x_b=y_b;
    end

    if len + norm(X-x1) > leng   % condition
        X = (X - x1) / norm(X-x1) * (leng - len) + x1;
        len = leng;
    else
        len = len + norm(X-x1);
    end
end
new_BC=cartesianToBarycentric(TR,new_ID,X);
end
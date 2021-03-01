function [x,y,x_b,y_b,vertex_x,vertex_y,edge_x,edge_y]=intersection(v1,v2,v3,a,b,c,d)
it1=[a b c]*v1';
it2=[a b c]*v2';
it3=[a b c]*v3';
ait1=abs(it1+d);
ait2=abs(it2+d);
ait3=abs(it3+d);
ep = 1e-15;
n_p=[a b c];
%vertex
if ait1<=ep && ait2>ep && ait3>ep && (it2+d)*(it3+d)>ep
    x=v1;x_b=0;
    vertex_x=1;
    y=v1;y_b=0;
    vertex_y=1;
    edge_x=0;edge_y=0;
elseif ait2<=ep && ait1>ep && ait3>ep && (it1+d)*(it3+d)>ep
    x=v2;x_b=0;
    vertex_x=2;
    y=v2;y_b=0;
    vertex_y=2;
    edge_x=0;edge_y=0;
elseif ait3<=ep && ait1>ep && ait2>ep && (it1+d)*(it2+d)>ep
    x=v3; x_b=0;
    vertex_x=3;
    y=v3; y_b=0;
    vertex_y=3;
    edge_x=0;edge_y=0;
%edge
elseif ait1<=ep && ait2<=ep && ait3>ep
    x=v1;x_b=0;
    vertex_x=1;
    y=v2;y_b=0;
    vertex_y=2;
    edge_x=0;edge_y=0;
elseif ait2<=ep && ait3<=ep && ait1>ep
    x=v2;x_b=0;
    vertex_x=2;
    y=v3;y_b=0;
    vertex_y=3;
    edge_x=0;edge_y=0;
elseif ait3<=ep && ait1<=ep && ait2>ep
    x=v1;x_b=0;
    vertex_x=1;
    y=v3;y_b=0;
    vertex_y=3;
    edge_x=0;edge_y=0;
%one vertex and one edge
elseif ait1<=ep && (it2+d)*(it3+d)<ep
    x=v1;x_b=0;
    vertex_x=1;
    s=-([a b c d]*[v2 1]')/(n_p*(v3-v2)');
    y=v2+s*(v3-v2);y_b=1;
    edge_y=[2,3];
    edge_x=0;vertex_y=0;
elseif ait2<=ep && (it1+d)*(it3+d)<ep
    x=v2;x_b=0;
    vertex_x=2;
    s=-([a b c d]*[v1 1]')/(n_p*(v3-v1)');
    y=v1+s*(v3-v1); y_b=1;
    edge_y=[1,3];
    edge_x=0;vertex_y=0;
elseif ait3<=ep && (it1+d)*(it2+d)<ep
    x=v3;x_b=0;
    vertex_x=3;
    s=-([a b c d]*[v1 1]')/(n_p*(v2-v1)');
    y=v1+s*(v2-v1); y_b=1; 
    edge_y=[1,2];
    edge_x=0;vertex_y=0;
% two edges
elseif (it1+d)*(it2+d)<ep && (it1+d)*(it3+d)<ep && (it2+d)*(it3+d)>ep
    s1=-([a b c d]*[v1 1]')/(n_p*(v2-v1)');
    x=v1+s1*(v2-v1);x_b=1;
    edge_x=[1,2];
    s2=-([a b c d]*[v1 1]')/(n_p*(v3-v1)');
    y=v1+s2*(v3-v1); y_b=1;
    edge_y=[1,3];
    vertex_x=0;vertex_y=0;
elseif (it1+d)*(it2+d)<ep && (it2+d)*(it3+d)<ep && (it1+d)*(it3+d)>ep
    s1=-([a b c d]*[v1 1]')/(n_p*(v2-v1)');
    x=v1+s1*(v2-v1);x_b=1;
    edge_x=[1,2];
    s2=-([a b c d]*[v2 1]')/(n_p*(v3-v2)');
    y=v2+s2*(v3-v2); y_b=1;
    edge_y=[2,3];
    vertex_x=0;vertex_y=0;
elseif (it1+d)*(it3+d)<ep && (it2+d)*(it3+d)<ep && (it1+d)*(it2+d)>ep
    s1=-([a b c d]*[v1 1]')/(n_p*(v3-v1)');
    x=v1+s1*(v3-v1);x_b=1;
    edge_x=[1,3];
    s2=-([a b c d]*[v2 1]')/(n_p*(v3-v2)');
    y=v2+s2*(v3-v2);y_b=1;
    edge_y=[2,3];
    vertex_x=0;vertex_y=0;
%face
%elseif ait1==0 && ait2==0 && ait3==0
end
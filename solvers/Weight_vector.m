function K = Weight_vector(w, n)

Weight_v=[];
a=[];
b=[];
for j=1:n-1
    weight=[];
    for i=1:n-j
        if w(j,i+j)~=0
            a=[a,j];
            b=[b,j+i];
            Weight_v=[Weight_v,w(j,i+j)];
        end
%         weight(1,i)=w(j,i+j);
    end
%     Weight_v=[Weight_v,weight];
end
G=graph(a,b);
k=-full(incidence(G));
K=k.*Weight_v;
K=K';

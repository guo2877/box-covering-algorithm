function [Boxes,Nb] = MEMB(Net,Rb)
% INPUTS: adjacency matrix, box radius
% OUTPUTS:box coverage, the number of boxes
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
tic;
D=Distance(Net);
Neighbor_Net=D;  
Neighbor_Net(Neighbor_Net>1)=0;    
D_Net=D;  
D_Net(D_Net>Rb)=0;  
D_Net(D_Net~=0)=1;   
Scale=size(Net,2);  
uncovered=1:Scale; 
covered=[];
centers=[];
while  ~isempty(uncovered)
    p=maximum_excluded_mass(covered,centers,D_Net);
    covered_p=find(D_Net(p,:)==1); 
    covered=union(covered,covered_p);
    covered=union(covered,p);
    centers=union(centers,p);
    uncovered=setdiff(uncovered,covered);
end

Nb=size(centers,2);
disp(num2str(Nb)); 
center_id = 1:Nb;
Boxes=zeros(1,Scale);
for i=1:Nb
   Boxes(1,centers(i))=center_id(i);
end
uncenters=setdiff(1:Scale,centers);
C_Net=D;   
C_Net(:,uncenters)=[];
dist=min(C_Net,[],2)';    
[~, index]=sort(dist);   
for m=1:Scale
    node=index(1,m);
    if dist(1,node)==0
        continue;
    else  
        neighbor=find(Neighbor_Net(node,:)==1);
        d_n=dist(neighbor);     
        min_nei=neighbor(d_n<dist(1,node));     
        box_min_nei=Boxes(1,min_nei);      
        Boxes(1,node)=box_min_nei(randperm(length(box_min_nei),1));    
    end
end
toc;
end

function [p] = maximum_excluded_mass(covered,centers,D_Net)
D_Net(:,covered)=0;
D_Net(centers,:)=0;
sum_excluded_mass=sum(D_Net,2);
max_excluded_mass=find(sum_excluded_mass==max(sum_excluded_mass));
p=max_excluded_mass(randperm(length(max_excluded_mass),1));
end
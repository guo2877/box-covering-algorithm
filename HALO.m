function [Boxes,Nb] = HALO(Net,Rb)
% INPUTS: adjacency matrix, box radius
% OUTPUTS:box coverage,the number of boxes
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
tic;
D=Distance(Net);
Scale=size(Net,2);  
Lb=2*Rb+1;
MC_cell=cell(1,1);  
runs=0;

D_Rb=D;
D_Rb(D_Rb>Rb)=0;
D_Rb(D_Rb~=0)=1; 
uncovered=1:Scale; 
covered=[];
centers=[];

while ~isempty(uncovered)
    p=maximum_excluded_mass(covered,centers,D_Rb);
    R=ball_of_seed_Rb(D_Rb,p);
    runs=runs+1;
    MC_cell{runs}=R; 
    covered=union(covered,R);
    centers=union(centers,p);
    uncovered=setdiff(uncovered,covered);
end

D_Lb=D;
D_Lb(D_Lb>=Lb)=0;
D_Lb(D_Lb~=0)=1;  
[deg,~,~]=degrees(Net);        
S=1:Scale;

while ~isempty(S)
    min_deg=find(deg==min(deg));
    seed=S(min_deg(randperm(length(min_deg),1))); 
    R=ball_of_seed_Lb(D_Lb,seed,S);
    for i=1:size(R,2)
        deg(S==R(i))=[];
        S(S==R(i))=[];
    end
    runs=runs+1;
    MC_cell{runs}=R;     
end
disp(num2str(runs)); 

MC=zeros(runs,Scale);
for i=1:runs
    MC(i,MC_cell{i})=1;   
end
DB=unique(MC,'rows','stable');  
freq=sum(DB,1);   
Boxes=zeros(1,size(DB,2));   
Color=1; 

r_sums=sum(DB,2);  
[~, index]=sort(r_sums);  
DB=DB(index,:);    
while ~isequal(freq,zeros(1,size(D,2)))
    cnode=find(DB(1,:)==1);
    fvalue=freq(cnode);  
    if size(find(fvalue==1),2)==0   
        freq(cnode)=freq(cnode)-1;  
        DB(1,:)=[];
    else
        Boxes(cnode)=Color;
        Color=Color+1;
        freq(cnode)=0;  
        DB(:,cnode)=0;
        DB(all(DB==0,2),:)=[];  
        r_sums=sum(DB,2);  
        [~, index]=sort(r_sums);   
        DB=DB(index,:);    
    end
end
Nb=max(Boxes);
toc;
end

function [R]=ball_of_seed_Lb(D,seed,S)
R=seed;
P=find(D(seed,:)==1);  
while ~isempty(P)
    pre_nodes=intersect(S,P);
    if ~isequal(size(pre_nodes,2),0)
        ui=pre_nodes(randperm(length(pre_nodes),1)); 
    else
        ui=P(randperm(length(P),1)); 
    end
    P=setdiff(P,ui);    
    R=union(R,ui);     
    P=intersect(P,find(D(ui,:)==1));  
end
end

function [p] = maximum_excluded_mass(covered,centers,D_Net)
D_Net(:,covered)=0;
D_Net(centers,:)=0;
sum_excluded_mass=sum(D_Net,2);
max_excluded_mass=find(sum_excluded_mass==max(sum_excluded_mass));   
p=max_excluded_mass(randperm(length(max_excluded_mass),1));
end

function [R]=ball_of_seed_Rb(D,seed)
R=seed;  
P=find(D(seed,:)==1); 
R=union(R,P); 
end






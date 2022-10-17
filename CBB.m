function [Boxes,Nb] = CBB(Net,Lb)
% INPUTS: adjacency matrix, box size
% OUTPUTS:box coverage, the number of boxes
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
tic;
D=Distance(Net);  
D(D>=Lb)=0;
D(D~=0)=1;   
Boxes=zeros(1,size(Net,2));   
Scale=size(Net,2);  
S=1:Scale;   
color=1;
while ~isempty(S)
    seed=S(randperm(length(S),1)); 
    box=ball_of_seed(D,seed,S);
    S=setdiff(S,box);
    Boxes(1,box)=color; 
    color=color+1;
end
Nb=max(Boxes);
toc;
end

function [box]=ball_of_seed(D,seed,S)
box=seed;
P=find(D(seed,:)==1); 
P=intersect(P,S);
while ~isempty(P)     
    ui=P(randperm(length(P),1));  
    P=setdiff(P,ui);
    box=union(box,ui);     
    P=intersect(P,find(D(ui,:)==1));  
end
end


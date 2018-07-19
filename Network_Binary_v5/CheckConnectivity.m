%% Connectivity Check

% Depth search
% visited = zero(1,913)
num_nodes = 652;
global visited Node
visited = zeros(1,num_nodes);
[Vessel, Node] = Create_Real_Rat_Mesentery();

dfs(20)

sum(visited)

function [] = dfs(n)
global visited Node   
    
    if visited(n)
        return
        
    else
        visited(n) = 1;
        for i = 1:Node{n}.num_Parent_Nodes
            if ~visited(Node{n}.Parent_Node(i))
                dfs(Node{n}.Parent_Node(i))
            end
                
        end
        for i = 1:Node{n}.num_Daughter_Nodes
            if ~visited(Node{n}.Daughter_Node(i))
                dfs(Node{n}.Daughter_Node(i))
            end
        end
    end
end


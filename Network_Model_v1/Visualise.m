function [] = Visualise(Vessel,Node)
%{
Function that recieves all the xy coordinates of the vessels/nodes and
plots it, for visualisation
%}

%{
Author = Michael Zhang
Date created = 11-06-18
v2
%}

global num_vessels num_nodes num_generations

for v = 1:num_vessels
    Coord_xy_Start(v,:) = Vessel{v}.xy_Start;
    Coord_xy_End(v,:) = Vessel{v}.xy_End;
    Length(v) = Vessel{v}.Length;
    Diameter(v) = Vessel{v}.Radius * 2;
    generation(v) = Vessel{v}.Generation;
end

x1 = Coord_xy_Start(:,1);
y1 = Coord_xy_Start(:,2);
x2 = Coord_xy_End(:,1);
y2 = Coord_xy_End(:,2);

%% Plot "Network"
figure
subplot(2,3,1)

plot([x1,x2]',[y1,y2]','--r')

%% Frequency vs Length

subplot(2,3,2)
histogram(Length,'Normalization','probability','BinWidth',50)
xlabel('Length')
ylabel('Frequency')
% Normalise maybe


%% Length vs Generation

length_avg = zeros(num_generations,1);
for g=1:num_generations
    length_avg(g) = mean(Length(generation == g));
    diameter_avg(g) = mean(Diameter(generation == g));
end

subplot(2,3,5)
% Averrage length data per generation
plot(1:num_generations,length_avg)
xlabel('Generation')
ylabel('Mean Length')

%% capillaries vs generation

subplot(2,3,4)

cap = zeros(num_nodes,1);
for n = 1:num_nodes
    if Node{n}.Terminal
        cap(n) = Node{n}.Generation;
    end
end

cap(cap == 0) = [];
histogram(cap)
xlabel('Generation')
ylabel('Number of Terminal segments')
% mean(cap)

subplot(2,3,3)
histogram(Diameter,'Normalization','probability','BinWidth',2.5)
xlabel('Diameter')
ylabel('Frequency')
% Normalise maybe

subplot(2,3,6)
plot(1:num_generations,diameter_avg)
xlabel('Generation')
ylabel('Mean Diameter')
end


function ID = getID(num_items, k)
    branches = 0:sum(num_items);
    for i = 1:size(num_items,2)
        ID{i} = branches((1+k:k+num_items(i)));
        k = k + num_items(i);
    end
end
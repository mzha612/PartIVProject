function [num] = Doubling(generations)
   % Could be improved
   num = 2*ones(1,2*generations);
   for i = 0:generations
       num(i+1)= num(i+1)^i;
   end
   num(generations+1:end) = fliplr(num(1:generations));
end
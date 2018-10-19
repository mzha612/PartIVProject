function [isValid] = Check_Valid2(Vessel,Node,sprout)


global num_vessels num_nodes
global x_min x_max y_min y_max


%
%
%
% for d = 1:2 % for each daughter vessel.
%
%
% end
L = 16.6;
for d = 1:2 % for each daughter vessel.
    %     check if intersecting any other exisiting vessel
    [p1,q1] = get_Coord(Vessel{Node{sprout}.Daughter_Vessel(d)});

    
    % extend length by a set amount
    %     ([p1.x q1.x], [p1.y q1.y])
    
    %     grad = (q1.y-p1.y)/(q1.x-p1.x);
%        norm([(q1.x-p1.x),(q1.y-p1.y)])
    theta = atan2d(q1.y-p1.y,q1.x-p1.x);
    q1.x = q1.x + L*cosd(theta);
    q1.y = q1.y + L*sind(theta);
%        norm([(q1.x-p1.x),(q1.y-p1.y)])
    %     (-1 + 2*(p1.x  < q1.x))
    %     (-1 + 2*(p1.y  < q1.y))
    %     check if intersecting any other exisiting vessel
%     [p1,q1] = get_Coord(Vessel{Node{sprout}.Daughter_Vessel(d)});
    if (p1.x < x_min) ||  (p1.x > x_max)|| (q1.x < x_min) ||  (q1.x > x_max)
        isValid = 0;
        return
    end
    
    
    if (p1.y < y_min) ||  (p1.y > y_max)|| (q1.y < y_min) ||  (q1.y > y_max)
        isValid = 0;
        return
    end
%     
    % Triangle
%     m = y_max/x_max;
%     
%     
%     
%     if (p1.y < y_min) ||  (p1.y > y_max)|| (q1.y < y_min) ||  (q1.y > y_max)
%         isValid = 0;
%         return
%     end
%     
%     if (p1.x < (p1.y-y_max)/m) ||  (p1.x > -(p1.y -y_max)/m)|| (q1.x < (q1.y-y_max)/m) ||  (q1.x > -(q1.y-y_max)/m)
%         isValid = 0;
%         return
%     end
%     
    for v = 1:num_vessels
        if (Node{sprout}.Parent_Vessel == v)
%                 || (Vessel{Node{sprout}.Daughter_Vessel(d)}.Daughter_Vessel(2) == v)
            continue;
        end
        [p2,q2] = get_Coord(Vessel{v});
        
        % extend length by a set amount
        %     ([p1.x q1.x], [p1.y q1.y])
        
        %     grad = (q1.y-p1.y)/(q1.x-p1.x);
%         norm([(q2.x-p2.x),(q2.y-p2.y)])
%          theta = atan2d(q2.y-p2.y,q2.x-p2.x);
%         q2.x = q2.x + L*cosd(theta);
%         q2.y = q2.y + L*sind(theta);
%         theta = atan2d(p2.y-q2.y,p2.x-q2.x);
%         p2.x = p2.x + L*cosd(theta);
%         p2.y = p2.y + L*sind(theta);
% %         norm([(q2.x-p2.x),(q2.y-p2.y)])
        if ((p1.x == p2.x) && (p1.y == p2.y)) || ...
                ((p1.x == q2.x) && (p1.y == q2.y)) || ((q1.x == q2.x) && (q1.y == q2.y))...
                || ((q1.x == p2.x) && (q1.y == p2.y))
            
           continue;
        end
        if (doIntersect(p1, q1, p2, q2))
            isValid = 0;
            return
        end
        
        
    end
    if (v == num_vessels) && (d == 2)
        isValid = 1;
        return
    end
    
end
% isValid = 0;
isValid = 1;

end


function [point1,point2] = get_Coord(Vessel)
point1.x = Vessel.xy_Start(1);
point1.y = Vessel.xy_Start(2);
point2.x = Vessel.xy_End(1);
point2.y = Vessel.xy_End(2);

end

function isOnSegment = onSegment(p,q,r)
if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) && ...
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
    isOnSegment = 1;
    return
end
isOnSegment = 0;
end


function ori = orientation(p,q,r)
val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

if (val == 0)
    ori = 0; % colinear
    return;
elseif (val > 0)
    ori = 1; % clock
    return;
else
    ori = 2; % counterclock
    return;
    
end

end


function isIntersecting = doIntersect(p1,q1,p2,q2)
o1 = orientation(p1, q1, p2);
o2 = orientation(p1, q1, q2);
o3 = orientation(p2, q2, p1);
o4 = orientation(p2, q2, q1);


%     // General case
if ((o1 ~= o2) && (o3 ~= o4))
    isIntersecting = 1;
    return;
end

%
%         // Special Cases
%         // p1, q1 and p2 are colinear and p2 lies on segment p1q1
if (o1 == 0 && onSegment(p1, p2, q1))
    isIntersecting = 1;
    return;
end

%             // p1, q1 and q2 are colinear and q2 lies on segment p1q1
if (o2 == 0 && onSegment(p1, q2, q1))
    isIntersecting = 1;
    return;
end
%                 // p2, q2 and p1 are colinear and p1 lies on segment p2q2
if (o3 == 0 && onSegment(p2, p1, q2))
    isIntersecting = 1;
    return;
end
%                    // p2, q2 and q1 are colinear and q1 lies on segment p2q2
if (o4 == 0 && onSegment(p2, q1, q2))
    isIntersecting = 1;
    return;
end

isIntersecting = 0;

%                         // Doesn't fall in any of the above cases

end
%
% function isIntersecting = doIntersect(p1,q1,p2,q2)
% distance(P1, Q1) = distance between P1 and Q1.
%
% dis =
%
% If point K lies on a line segment then following equation is true : distance(P1,Q1) = distance(P1,K) + distance(K,Q1)
%
% Find the slopes of two line segments and let's say if they are not equal - then
%
% 1) Form two equations of the form ax+by+c = 0
% 2) Compute the intersecting point by solving above line equations.Let's call our intersecting point as I.
% 3) Verify whether distance(P1, Q1) == distance(P1,I) + distance(I,Q1)
% 4) Verify the same with P2 and Q2 as well
% 5) If intersecting point lies on both the line segments then above two verifications should succeed
%
% If slopes are equal then check if their equations are same or not. If equations are different then they are parallel lines - no intersecting point
%
% If equations are same then the points are collinear - verify the following (one of the points of the line segment should lie on the other)
%
% distance(P1,Q1) == distance(P1,P2) + distance(P2, Q1)
% distance(P1,Q1) == distance(P1,Q2) + distance(Q2, Q1)
%
% If any of the above statements are true then the line segments intersect
% end
function check_max_principle(u_aux,uFCT)
% checks max principle at each node

% size of system
n = length(u_aux);

% check max principle at each node, except for Dirichlet node, where
% technically, the max principle won't necessarily be satisfied since
% R+ and R- are forced to 1
for i = 2:n
    % compute index range of support of i
    i1=max(i-1,1);
    i2=min(i+1,n);
    % compute max and min old solution in the support of each dof
    u_max = max(u_aux(i1:i2));
    u_min = min(u_aux(i1:i2));
    % check max principle
    if (uFCT(i) > (u_max+eps))
        error('Did not satisfy max principle at i = %i\n',i);
    end
    if (uFCT(i) < (u_min-eps))
        error('Did not satisfy max principle at i = %i\n',i);
    end
end

return
end
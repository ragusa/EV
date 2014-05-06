function check_max_principle(u_new,W_max,W_min)
% checks max principle at each node
%
% u_new = new solution
% W_max = upper bounds for max principle
% W_min = upper bounds for max principle

% size of system
n = length(u_new);

% check max principle at each node, except for Dirichlet node, where
% technically, the max principle won't necessarily be satisfied since
% R+ and R- are forced to 1
for i = 2:n
    if (u_new(i) > (W_max(i)+eps))
        warning('Exceeded upper bound of max principle at node %i\n',i);
    end
    if (u_new(i) < (W_min(i)-eps))
        warning('Exceeded upper bound of max principle at node %i\n',i);
    end
end

return
end
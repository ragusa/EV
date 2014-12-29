function check_DMP(u,Wplus,Wminus,description)

n = length(u);

for i = 1:n    
    if (u(i) > Wplus(i))
        fprintf('%s:\tViolated upper bound for node %i: W+ = %f, u = %f\n',...
            description,i,Wplus(i),u(i));
    elseif (u(i) < Wminus(i))
        fprintf('%s:\tViolated lower bound for node %i: W- = %f, u = %f\n',...
            description,i,Wminus(i),u(i));
    end
end

end
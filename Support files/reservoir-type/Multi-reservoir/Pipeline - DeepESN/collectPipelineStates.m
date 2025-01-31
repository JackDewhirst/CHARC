function[final_states,individual] = collectPipelineStates(individual,input_sequence,config)    

%% Collect states for plain ESN
%if single input entry, add previous state
if size(input_sequence,1) == 1
    input_sequence = [zeros(size(input_sequence)); input_sequence];
end

for i= 1:config.num_reservoirs
    if size(input_sequence,1) == 2
        states{i} = individual.last_state{i};
    else
        states{i} = zeros(size(input_sequence,1),individual.nodes(i));
    end
    x{i} = zeros(size(input_sequence,1),individual.nodes(i));
end


%equation: x(n) = f(Win*u(n) + S)
for n = 2:size(input_sequence,1)
    
    for i= 1:config.num_reservoirs
        
        x{i}(n,:) = x{i}(n,:) + ((individual.W{i,i}*individual.W_scaling(i,i))'*states{i}(n-1,:)')';

        if i == 1
            input = ([individual.bias_node input_sequence(n,:)])';
        else
            input = ([individual.bias_node input_sequence(n,:) states{i-1}(n,:)])';
        end
        
        if iscell(individual.activ_Fcn)
            for p = 1:individual.nodes(i)            
                states{i}(n,p) = feval(individual.activ_Fcn{p},((individual.input_weights{i}(p,:)*individual.input_scaling(i))*input) + x{i}(n,p)'); 
            end
        else
            states{i}(n,:) = feval(individual.activ_Fcn,((individual.input_weights{i}*individual.input_scaling(i))*input) + x{i}(n,:)'); 
        end
        
    end
    
end


if config.leak_on
    states = getLeakStates(states,individual,input_sequence,config);
end

% concat all states for output weights
final_states = [];
for i= 1:config.num_reservoirs
    final_states = [final_states states{i}];
    
    %assign last state variable
    individual.last_state{i} = states{i}(end,:);
end

if config.add_input_states == 1
    final_states = [final_states input_sequence];
end

final_states = final_states(config.wash_out+1:end,:); % remove washout
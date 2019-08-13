%% Mutation operator used for different reservoir systems
% Details:
% - number of weights mutated is based on mut_rate; 50% chance to change existing weight or remove it
% - 25% chance to change global parameters


%-
%J.D; Important note
%if config.connectiondensity = 0 and there are sub-reservoirs
%then this function will break
%-

function offspring = mutateRoR(offspring,config)

% params - input scaling and leak rate
input_scaling = offspring.input_scaling(:);
pos =  randi([1 length(input_scaling)],rand < config.mut_rate,1);
input_scaling(pos) = 2*rand(length(pos),1)-1;
offspring.input_scaling = reshape(input_scaling,size(offspring.input_scaling));

leak_rate = offspring.leak_rate(:);
pos =  randi([1 length(leak_rate)],rand < config.mut_rate,1); % 0.25
leak_rate(pos) = rand(length(pos),1);
offspring.input_scaling = reshape(leak_rate,size(offspring.leak_rate));

% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs

    % input weights
    input_weights = offspring.input_weights{i};
    pos =  randi([1 length(input_weights)],ceil(config.mut_rate*length(input_weights)),1);
    for n = 1:length(pos)
        if rand < 0.5 % 50% chance to zero weight
            input_weights(pos(n)) = 0;
        else
            input_weights(pos(n)) = 2*rand-1;
        end
    end
    offspring.input_weights{i} = reshape(input_weights,size(offspring.input_weights{i}));
        
    % W scaling
    W_scaling = offspring.W_scaling(i,:);
    pos =  randi([1 length(W_scaling)],rand < config.mut_rate,1);
    W_scaling(pos) = 2*rand(length(pos),1);
    offspring.W_scaling(i,:) = reshape(W_scaling,size(offspring.W_scaling(i,:)));

    % hidden weights
    for j = 1:config.num_reservoirs
        W = offspring.W{i,j}(:);
        
        if i == j
             % select weights to change
            pos =  randi([1 length(W)],ceil(config.mut_rate*length(W)),1);
            for n = 1:length(pos)
                if rand < 0.5 % 50% chance to zero weight
                   W(pos(n)) = 0;
                else
                    W(pos(n)) = rand-0.5;
                end   
            end
                
        else %the result of the following is that the size of the connection matrix ... 
            %doesn't change during mutation, hopefully
            pos =  randi([1 length(W)],ceil(config.mut_rate*length(W)),1);
            for n = 1:length(pos)
                non_zero = find(W);
                if ismember(pos(n),non_zero) %if position is non_zero
                    if rand < 0.5 %50 chance to make position zero
                        W(pos(n)) = 0;
                        not_found = 1;
                        while(not_found)  
                            flip = randi([1 length(W)]); %find another position
                            if ~ismember(flip,non_zero) %check if it's zero
                                not_found = 0;
                            end
                        end
                        W(flip) = rand-0.5; % if it is, turn it on
                    end
                    
                else                        %if position is zero
                    if rand < 0.5 %50 chance to make position non_zero
                        non_zero = find(W);
                        W(pos(n)) = rand-0.5;
                        %not_found = 1;
                        flip = randi([1 size(non_zero,2)]); %find another position
                        non_zero = find(W);
                        flip = non_zero(flip);
                        %while(not_found)  
                        %    flip = randi([1 length(nonzero)]); %find another position
                        %    flip = nonzero(flip);
                        %    if ismember(flip,non_zero) %check if it's non_zero
                        %        not_found = 0;
                        %    end
                        %end
                        W(flip) = 0; % if it is, turn it off
                    end
                end
            end
            
            %
            if( nnz(W) < ceil(config.connection_density ...
                    *(offspring.nodes(i)^2 )) )
                
                %pick a zero element, add more elements until we hit target
                while( nnz(W) < ceil(config.connection_density ...
                    *(offspring.nodes(i)^2 )) )
                    W(randperm(length(W),1)) = rand - 0.5; 
                end
                
            elseif( nnz(W) > ceil(config.connection_density ...
                    *(offspring.nodes(i)^2 )) )
                
                %pick a non-zero element, delete it until we target
                while( nnz(W) > ceil(config.connection_density ...
                    *(offspring.nodes(i)^2 )) )
                
                    W(randperm(length(W),1)) = 0;
                    
                end
                
            else
                
            end
                
            
            %while (density < config.connection_density - 0.05)  || ...
            %         (density > config.connection_density + 0.05)
            %    pos =  randi([1 length(W)],ceil(config.mut_rate*length(W)),1);
            %    for n = 1:length(pos)
            %        if rand < 0.5 % 50% chance to zero weight
            %            W(pos(n)) = 0;
            %        else
            %            W(pos(n)) = rand-0.5;
            %        end   
            %    end
            %    density = nnz(W)/ (size(W,2)*size(W,1));
            %end
            
            offspring.W{i,j} = reshape(W,size(offspring.W{i,j}));
            
        end
            
        
    end
end

% mutate output weights
if config.evolve_output_weights
    output_weights = offspring.output_weights(:);
    pos =  randi([1 length(output_weights)],ceil(config.mut_rate*length(output_weights)),1);
     for n = 1:length(pos)
        if rand > 0.75 % 75% chance to zero weight
            output_weights(pos(n)) = 0;
        else
            output_weights(pos(n)) = 2*rand-1;
        end
    end
    offspring.output_weights = reshape(output_weights,size(offspring.output_weights));
end



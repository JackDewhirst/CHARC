function population = createRoR(config)


%% Reservoir Parameters
for pop_indx = 1:config.pop_size

    % add performance records
    population(pop_indx).train_error = 1;
    population(pop_indx).val_error = 1;
    population(pop_indx).test_error = 1;
    
    % add single bias node
    population(pop_indx).bias_node = 1;
    
    % assign input/output count
    if isempty(config.train_input_sequence) 
        population(pop_indx).n_input_units = 1;
        population(pop_indx).n_output_units = 1;
    else
        population(pop_indx).n_input_units = size(config.train_input_sequence,2);
        population(pop_indx).n_output_units = size(config.train_output_sequence,2);
    end
    
    % iterate through subreservoirs
    for i = 1:config.num_reservoirs
        
        %define num of units
        population(pop_indx).nodes(i) = config.num_nodes(i);

        
        % Scaling and leak rate
        population(pop_indx).input_scaling(i) = 2*rand-1; %increases nonlinearity
        population(pop_indx).leak_rate(i) = rand;
        
       %inputweights
        if config.sparse_input_weights
            input_weights = sprand(population(pop_indx).nodes(i),  population(pop_indx).n_input_units+1, 0.1); 
            input_weights(input_weights ~= 0) = ...
                2*input_weights(input_weights ~= 0)  - 1;
            population(pop_indx).input_weights{i} = input_weights;
        else
            population(pop_indx).input_weights{i} = 2*rand(population(pop_indx).nodes(i),  population(pop_indx).n_input_units+1)-1; 
        end
        
              
        %assign different activations, if necessary
        if config.multi_activ 
            activ_positions = randi(length(config.activ_list),1,population(pop_indx).nodes(i));
            for act = 1:length(activ_positions)
                population(pop_indx).activ_Fcn{i,act} = config.activ_list{activ_positions(act)};
            end
        else
            population(pop_indx).activ_Fcn = config.activ_list;
        end
        
        population(pop_indx).last_state{i} = zeros(1,population(pop_indx).nodes(i));
    end
    
    %track total nodes in use
    population(pop_indx).total_units = 0;
    
    
    %% weights and connectivity of all reservoirs
    % I want to set the density to some fixed valued for the connection matrices,
    % then I also want to go into mutateRoR and force the connection
    % matrices to stay at the same value of denisty throughout the process
    for i= 1:config.num_reservoirs
        
        for j= 1:config.num_reservoirs 
            if config.connection_density == 0 || i==j
                population(pop_indx).connectivity(i,j) =  10/population(pop_indx).nodes(i); 

                internal_weights = sprand(population(pop_indx).nodes(i), population(pop_indx).nodes(j), population(pop_indx).connectivity(i,j));
                internal_weights(internal_weights ~= 0) = ...
                    internal_weights(internal_weights ~= 0)  - 0.5;

                % assign scaling for inner weights 
                population(pop_indx).W_scaling(i,j) = 2*rand;            
                population(pop_indx).W{i,j} = internal_weights; 
            else %I want for the below to just take config.connection_density and give me that many conections
                population(pop_indx).connectivity(i,j) =  config.connection_density;
                pos = randperm(population(pop_indx).nodes(i)^2, ceil(config.connection_density ...
                    *(population(pop_indx).nodes(i))^2));
                internal_weights = sparse(population(pop_indx).nodes(i),population(pop_indx).nodes(i));
                %internal_weights = sparse(population(pop_indx).nodes(i),population(pop_indx).nodes(i));
                for k=1:length(pos)
                    internal_weights(pos(k)) = rand - 0.5;
                end

                % assign scaling for inner weights 
                population(pop_indx).W_scaling(i,j) = 2*rand;            
                population(pop_indx).W{i,j} = internal_weights; 
            end

        end
        
        %I then want to prune the interconnects. e.g. two 50 node
        %reservoirs are gonna be in a 100 node matrix. areas
        %(50-100,1-50) and (1-50,50-100) represent connections between the
        %two reservoirs. I can fiddle the 
        
        population(pop_indx).total_units = population(pop_indx).total_units + population(pop_indx).nodes(i); 
    end
    

    % add rand output weights
    if config.add_input_states
        population(pop_indx).output_weights = 2*rand(population(pop_indx).total_units + population(pop_indx).n_input_units, population(pop_indx).n_output_units)-1;      
    else
        population(pop_indx).output_weights = 2*rand(population(pop_indx).total_units, population(pop_indx).n_output_units)-1;
    end
    
    population(pop_indx).behaviours = [];
    
end
%% Select Data Script: Generate task data sets and split data
function [config] = selectDataset(config)

scurr = rng;
temp_seed = scurr.Seed;

rng(1,'twister');

switch config.dataset
    
    %% Chaotic systems
    case 'NARMA10' %input error 4 - good task
        err_type = 'NMSE';
        queue_type = 'simple';
        wash_out =200;
        sequence_length = 8000;
        train_fraction=0.25;    val_fraction=0.375;    test_fraction=0.375;
        [input_sequence,output_sequence] = generate_new_NARMA_sequence(sequence_length,10);
        input_sequence = 2*input_sequence-0.5;
        output_sequence = 2*output_sequence-0.5;
        %fprintf('NARMA 10 task: %s \n',datestr(now, 'HH:MM:SS'))
        
    case 'NARMA20' %input error 4 - good task
        err_type = 'NMSE';
        queue_type = 'simple';
        wash_out =200;
        sequence_length = 8000;
        train_fraction=0.25;    val_fraction=0.375;    test_fraction=0.375;
        [input_sequence,output_sequence] = generate_new_NARMA_sequence(sequence_length,20);
        input_sequence = 2*input_sequence-0.5;
        output_sequence = 2*output_sequence-0.5;
        %fprintf('NARMA 20 task: %s \n',datestr(now, 'HH:MM:SS'))
        
    case 'NARMA30' %input error 4 - good task
        err_type = 'NMSE';
        queue_type = 'simple';
        wash_out =200;
        sequence_length = 8000;
        train_fraction=0.25;    val_fraction=0.375;    test_fraction=0.375;
        [input_sequence,output_sequence] = generate_new_NARMA_sequence(sequence_length,30);
        input_sequence = 2*input_sequence-0.5;
        output_sequence = 2*output_sequence-0.5;
        %fprintf('NARMA 30 task: %s \n',datestr(now, 'HH:MM:SS'))
        
    case 'NARMA10_DLexample' %input error 4 - good task
        err_type = 'NRMSE';
        queue_type = 'simple';
        config.preprocess = 0;
        
        wash_out = 100;
        sequence_length = 2700;
        train_fraction = 0.3333;    val_fraction=0.3333;    test_fraction=0.3333;
        
        [input_sequence,output_sequence] = generate_new_NARMA_sequence(sequence_length,10);
        
    case 'NARMA10_QRC' %input error 4 - good task
        err_type = 'NSE';
        queue_type = 'simple';
        config.preprocess = 0;
        
        wash_out = 2000;
        sequence_length = 12000;
        train_fraction = 0.3333;    val_fraction=0.3333;    test_fraction=0.3333;
        
        [input_sequence,output_sequence] = generate_new_NARMA_sequence(sequence_length,10);
        input_sequence = (input_sequence*2)*0.2;
        
    case 'HenonMap' % input error > 1 - good task
        queue_type = 'simple';
        err_type = 'NMSE';
        wash_out =200;
        sequence_length= 8000;
        stdev = 0.05;
        train_fraction=0.25;    val_fraction=0.375;    test_fraction=0.375;
        [input_sequence,output_sequence] = generateHenonMap(sequence_length,stdev);
        
        %% Time-series
    case 'IPIX_plus5' % good task
        err_type = 'IPIX';
        queue_type = 'Weighted';
        wash_out =100;
        sequence_length = 2000;
        train_fraction=0.4;    val_fraction=0.25;    test_fraction=0.35;   %val and test are switched later so ratios need to be swapped
        
        % IPIX radar task
        %load hiIPIX.txt
        load loIPIX.txt
        
        ahead = 5;
        data = loIPIX(1:sequence_length+ahead,:);
        input_sequence = data(1:sequence_length,:);
        output_sequence = data(ahead+1:end,:);
        
        %fprintf('Low IPIX task - 5 ahead. \n Started at %s \n',datestr(now, 'HH:MM:SS'))
        
    case 'IPIX_plus1' % good task
        err_type = 'IPIX';
        queue_type = 'Weighted';
        wash_out =100;
        sequence_length = 2000;
        train_fraction=0.4;    val_fraction=0.25;    test_fraction=0.35;   %val and test are switched later so ratios need to be swapped
        
        % IPIX radar task
        %load hiIPIX.txt
        load loIPIX.txt
        
        ahead = 1;
        data = loIPIX(1:sequence_length+ahead,:);
        input_sequence = data(1:sequence_length,:);
        output_sequence = data(ahead+1:end,:);
        
        %fprintf('Low IPIX task 1 ahead. \n Started at %s \n',datestr(now, 'HH:MM:SS'))
        
        
    case 'Laser' % good task
        queue_type = 'simple';
        err_type = 'NMSE';
        % Sante Fe Laser generator task
        wash_out =200;
        sequence_length = 8000;
        train_fraction=0.25;    val_fraction=0.375;    test_fraction=0.375;
        
        ahead = 1;
        data = laser_dataset;  %checkout the list at http://uk.mathworks.com/help/nnet/gs/neural-network-toolbox-sample-data-sets.html
        data = cell2mat(data(:,1:sequence_length+ahead));
        input_sequence = data(1:end-ahead)';
        output_sequence = data(ahead+1:end)';
        
        %fprintf('Laser task TSP - 64 electrode test: %s \n',datestr(now, 'HH:MM:SS'))
        
        
    case 'Sunspot' % good task but not sure about dataset- problem with dividing set
        queue_type = 'simple';
        err_type = 'NMSE';
        % Sunspot task - needs proper dataset separation
        wash_out =100;
        sequence_length = 3100;
        train_fraction= 1600/sequence_length;    val_fraction=500/sequence_length;    test_fraction=1000/sequence_length;
        
        ahead = 1;
        load sunspot.txt %solar_dataset;  %checkout the list at http://uk.mathworks.com/help/nnet/gs/neural-network-toolbox-sample-data-sets.html
        data = sunspot(1:sequence_length+ahead,4);
        input_sequence = data(1:end-ahead);
        output_sequence = data(ahead+1:end);
        
        %fprintf('Sunspot task TSP: %s \n',datestr(now, 'HH:MM:SS'))
        
        %% Pattern Recognition - using PCA to reduce dimensions maybe very useful
        
    case 'Autoencoder'
        err_type = 'MSE';
        queue_type = 'simple';
        wash_out = 0;
        train_fraction=0.25;    val_fraction=0.375;    test_fraction=0.375;
        
        t = digitTrainCellArrayData; %28 x 28 image x 5000
        for i=1:length(t)
            u(:,i) = t{i}(:);
        end
        
        input_sequence= u';
        output_sequence= u';
        
    case 'NIST-64' %Paper: Reservoir-based techniques for speech recognition
        err_type = 'OneVsAll_NIST';
        queue_type = 'Weighted';
        wash_out =150;
        xvalDetails.kfold = 5;
        xvalDetails.kfoldSize = 150;
        xvalDetails.kfoldType = 'standard';
        train_fraction=0.7;    val_fraction=0.15;    test_fraction=0.15; %meaningless
        
        y_list = [];
        u_list = [];
        lens = [];
        
        for i = 1:5
            l = [ 1 2 5 6 7];
            for j = 1:10
                for n = 0:9
                    u_z = zeros(77,xvalDetails.kfoldSize);
                    u = load(strcat('s',num2str(l(i)),'_u',num2str(j),'_d',num2str(n)));
                    u_z(:,1:size(u.spec,2)) = u.spec;
                    
                    y = zeros(10,size(u_z,2))-1;
                    y(n+1,:) = ones(1,size(u_z,2));
                    u_list = [u_list u_z];
                    lens = [lens size(u_z,2)];
                    y_list = [y_list y];
                end
            end
        end
        
        input_sequence = u_list';
        output_sequence = y_list';
        
        
    case 'NonChanEqRodan' % (1:in, 1:out) error 0.999 Good task, requires memory
        err_type = 'NMSE';
        queue_type = 'simple'; %input alone error = 0.091
        wash_out =200;
        sequence_length = 8000;
        train_fraction=0.25;    val_fraction=0.375;    test_fraction=0.375;
        
        [input_sequence, output_sequence] = NonLinear_ChanEQ_data(sequence_length);
        input_sequence =input_sequence';
        output_sequence =output_sequence';
        
    case 'handDigits'
        
        err_type = 'softmax';
        queue_type = 'Weighted';
        wash_out =10;
        train_fraction=0.8;    val_fraction=0.1;    test_fraction=0.1;
        dataset_length = 5000; %manually change dataset length for xval
        
        load('handDigits.mat');
        input_sequence = X;
        output_sequence = [];
        for i = 1:10
            output_sequence(:,i) = y==i;
        end
        
        target=randperm(dataset_length);
        temp_inputSequence = input_sequence(target,:);
        temp_outputSequence = output_sequence(target,:);
        
        input_sequence = temp_inputSequence;
        output_sequence = temp_outputSequence;
        
    case 'JapVowels' %(12: IN, 9:OUT - binary ) - input only 83% accuracy!  Train:0.2288  Test:0.1863
        err_type = 'softmax'; %Paper: Optimization and applications of echo state networks with leaky- integrator neurons
        queue_type = 'Weighted';
        % Nine male speakers uttered two Japanese vowels /ae/ successively.
        % For each utterance, with the analysis parameters described below, we applied
        % 12-degree linear prediction analysis to it to obtain a discrete-time series
        % with 12 LPC cepstrum coefficients. This means that one utterance by a speaker
        % forms a time series whose length is in the range 7-29 and each point of a time
        % series is of 12 features (12 coefficients).
        % The number of the time series is 640 in total. We used one set of 270 time series for
        % training and the other set of 370 time series for testing.
        wash_out =100;
        
        [train_input_sequence,trainOutputSequence,testInputSequence,test_output_sequence] = readJapVowels();
        input_sequence = [train_input_sequence; testInputSequence];
        output_sequence = [trainOutputSequence; test_output_sequence];
        train_fraction=size(train_input_sequence,1)/9961;    val_fraction=(size(testInputSequence,1)/9961)*0.1;    test_fraction=(size(testInputSequence,1)/9961)*0.9;
        
        t =  randperm(dataset_length,dataset_length);
        
    case 'SignalClassification'
        err_type = 'softmax';
        queue_type = 'simple';
        wash_out =100;
        train_fraction=0.25;    val_fraction=0.375;    test_fraction=0.375;
        
        freq = 1000;
        fprintf('Signal Classification: \n',datestr(now, 'HH:MM:SS'))
        fprintf('Freq: %d Hz\n',freq);
        scanFreq = 20000; %per channel
        step = 1/scanFreq;
        t = 0:step:1-step;
        amplitude = 1;
        sequence_length = 6000;
        period = 20;
        
        % sinewave input
        input_sequence(:,1) = amplitude*sin(2*pi*freq*t);
        input_sequence(:,2) = amplitude*square(2*pi*freq*t);
        
        cnt = 1; sinInput =[];squareInput=[];
        for i = 0:period:sequence_length-period
            sinInput(cnt,i+1:i+period) = input_sequence(i+1:i+period,1);
            squareInput(cnt,i+1:i+period) = input_sequence(i+1:i+period,2);
            cnt = cnt +1;
        end
        
        combInput = zeros(sequence_length,1);
        combOutput= ones(sequence_length,2)*0;
        for i = 1:sequence_length/period
            if round(rand)
                combInput = combInput+sinInput(i,:)';
                combOutput((i*period)-period+1:i*period,1) =  ones(period,1);
            else
                combInput = combInput+squareInput(i,:)';
                combOutput((i*period)-period+1:i*period,2) =  ones(period,1);
            end
        end
        
        input_sequence = combInput;
        output_sequence = combOutput;
        
    case 'Iris' %iris_dataset; (4:in, 3:out) %input alone 76% - medium task
        err_type = 'IJCNNpaper';%'IJCNNpaper';%'confusion';
        queue_type = 'Weighted';
        wash_out = 0;
        % train_fraction=0.66666667;    val_fraction=0.333333/2;    test_fraction=0.333333/2;
        train_fraction=0.5;    val_fraction=0.25;    test_fraction=0.25;
        dataset_length = 150;
        
        t =  randperm(dataset_length,dataset_length);
        
        [input_sequence, output_sequence] =  iris_dataset; %iris_dataset; (4:in, 3:out)
        input_sequence = input_sequence(:,t)';
        output_sequence = output_sequence(:,t)';
        
    case 'MSO'
        err_type = 'NRMSE';
        queue_type = 'simple'; %?
        wash_out =100;
        sequence_length= 2000;
        train_fraction=0.583333;    val_fraction=0.16667;    test_fraction=0.25;
        
        for t = 1:sequence_length
            u(1,t) = sin(0.2*t)+sin(0.311*t);
            u(2,t) = sin(0.2*t)+sin(0.311*t)+sin(0.42*t);
            u(3,t) = sin(0.2*t)+sin(0.311*t)+sin(0.42*t)+sin(0.51*t);
            u(4,t) = sin(0.2*t)+sin(0.311*t)+sin(0.42*t)+sin(0.51*t)+sin(0.74*t);
            %u(1,t) = sin(0.2*t)+sin(0.311*t)+sin(0.42*t)+sin(0.51*t)+sin(0.63*t)+sin(0.74*t)+sin(0.85*t)+sin(0.97*t);
        end
        %predictor - not sure what predictor value is best
        %ahead = 10;
        %outputSequence = u(:,ahead+1:end)';
        input_sequence = zeros(sequence_length,size(u,1));%u(:,1:end-ahead)';
        output_sequence = u';
        
    case 'secondOrderTask' %best 3.61e-3
        queue_type = 'simple';
        err_type = 'NMSE';
        
        wash_out =50;
        sequence_length = 700;
        train_fraction= 300/sequence_length;    val_fraction=100/sequence_length;    test_fraction=300/sequence_length;
        u = rand(sequence_length,1)/2;
        y = zeros(sequence_length,1);
        for i = 3:sequence_length
            y(i) = 0.4*y(i-1)+0.4*y(i-1)*y(i-2)+0.6*(u(i).^2) + 0.1;
        end
        input_sequence = u;
        output_sequence = y;
        
    case 'MNIST'
        err_type = 'softmax';
        queue_type = 'simple';
        wash_out = 0;
        train_fraction=0.25;    val_fraction=0.375;    test_fraction=0.375;
        preprocess = 0;
        
        [xTrainImages,tTrain] = digitTrainCellArrayData; %28 x 28 image x 5000
        
        for i=1:length(xTrainImages)
            u(:,i) = xTrainImages{i}(:);
        end
        
        target=randperm(length(xTrainImages));
        temp_inputSequence = u(:,target);
        temp_outputSequence = tTrain(:,target);
        
        input_sequence = temp_inputSequence';
        output_sequence = temp_outputSequence';
        
        
    case 'poleBalance'
        err_type = 'empty';
        queue_type = 'empty';
        train_fraction=0.1;    val_fraction=0.1;    test_fraction=0.1;
        wash_out = 1;
        input_sequence= zeros(100,6);
        output_sequence= zeros(100,1);
        
        
    case 'robot'
        err_type = 'empty';
        queue_type = 'empty';
        train_fraction=0.1;    val_fraction=0.1;    test_fraction=0.1;
        wash_out = 1;
        config.num_sensors = 8;
        input_sequence= zeros(100,config.num_sensors+1);
        output_sequence= zeros(100,4);
        
    case 'CPPN'
        err_type = 'empty';
        queue_type = 'empty';
        train_fraction=0.1;    val_fraction=0.1;    test_fraction=0.1;
        wash_out = 1;
        input_sequence= zeros(100,config.CPPN_inputs);
        output_sequence= zeros(100,config.CPPN_outputs);
               
    case 'BinaryNbitAdder'
        
        err_type = 'hamming';
        queue_type = 'simple';
        type = 'nbit_adder';
        bit = 3;
        datalength = 5000;
        wash_out = 25;
        train_fraction=0.5;    val_fraction=0.25;    test_fraction=0.25;
        config.preprocess =0;
        
        A_in = randi([0 (2^bit)-1],datalength,1);
        B_in = randi([0 (2^bit)-1],datalength,1);
        
        input = [de2bi(A_in,bit) de2bi(B_in,bit)];
        output=[];
        for i = 1:datalength
            output(i,:) = getAdderTruthTable(type,[bit,A_in(i),B_in(i)]);
        end
        
        in  = [bi2de(input(:,1:bit)) bi2de(input(:,bit+1:end))];
        out = bi2de(output);
        hist(out)
        
        % get uniform distribution
        [N,edges,bin] = histcounts(out,2^bit*2);
        cnt = 1;bin_in=[];bin_out=[];
        for i = min(bin):max(bin)
            bin_in{cnt} = input(bin == i,:);
            bin_out{cnt} = output(bin == i,:);
            cnt = cnt +1;
        end
        
        for i = 1:datalength
            ex = 1;
            while(ex)
                pos = randi([min(bin) max(bin)]);
                if ~isempty(bin_in{pos})
                    ex = 0;
                end
            end
            pos2 = randi([1 length(bin_in{pos})]);
            input_sequence(i,:) = bin_in{pos}(pos2,:);
            output_sequence(i,:) = bin_out{pos}(pos2,:);
        end
        
        in  = [bi2de(input_sequence(:,1:bit)) bi2de(input_sequence(:,bit+1:end))];
        out = bi2de(output_sequence);
        hist(out)
        %hist(in(:,2))
        %hist(in(:,1))
        
    case 'ImageGaussian' % Gaussian noise task
        err_type = 'softmax';
        queue_type = 'simple';
        
        wash_out = 0;
        train_fraction=0.5;    val_fraction=0.25;    test_fraction=0.25;
        config.preprocess = 0;
        
        startscript;
        
        imagesCombined = horzcat(imagesOriginal,imagesGaussian);
        
        inputs =[];trainingTarget=[];
        for i=1:length(imagesCombined)
            %inputs = imagesCombined{i}];
            
            if (i <= (length(imagesCombined)/2))
                %trainingTarget(:,i,1) = 0;
                %trainingTarget(:,i,2) = 1;
                trainingTarget{i} = repmat([0 1],size(imagesCombined{i},1),1);
                
            else
                %trainingTarget(:,i,1) = 1;
                %trainingTarget(:,i,2) = 0;
                trainingTarget{i} = repmat([1 0],size(imagesCombined{i},1),1);
            end
            
        end
        
        target=randperm(length(imagesCombined));
        temp_inputSequence = [];temp_outputSequence=[];
        for i = 1:length(target)
            temp_inputSequence = [temp_inputSequence; imagesCombined{target(i)}];
            temp_outputSequence = [temp_outputSequence; trainingTarget{target(i)}];
        end
        
        input_sequence = temp_inputSequence;
        output_sequence = temp_outputSequence;
end

%% preprocessing
if config.preprocess
    for i = 1:size(input_sequence,2)
        input_sequence(input_sequence(:,i) ~= 0,i) = (input_sequence(input_sequence(:,i) ~= 0,i)-mean(input_sequence(:,i)))/(max(input_sequence(:,i))-min(input_sequence(:,i)));
    end
    
    for i = 1:size(output_sequence,2)
        output_sequence(output_sequence(:,i) ~= 0,i) = (output_sequence(output_sequence(:,i) ~= 0,i)-mean(output_sequence(:,i)))/(max(output_sequence(:,i))-min(output_sequence(:,i)));
    end
end

if config.discrete %choose n-bit word length if needed by adding s,w,f to func() parameters
    
    if config.parallel
        config.poolobj = gcp;
        addAttachedFiles(config.poolobj,{'bin2num.m'})
    end
    
    [input_sequence, config.q] = double2binaryInputVector(input_sequence,config.nbits);
    [output_sequence, config.q] = double2binaryInputVector(output_sequence,config.nbits);
end

[train_input_sequence,val_input_sequence,test_input_sequence] = ...
    split_train_test3way(input_sequence,train_fraction,val_fraction,test_fraction);

[train_output_sequence,val_output_sequence,test_output_sequence] = ...
    split_train_test3way(output_sequence,train_fraction,val_fraction,test_fraction);

% Go back to old seed
rng(temp_seed,'twister');

% squash into structure
config.train_input_sequence = train_input_sequence;
config.train_output_sequence = train_output_sequence;
config.val_input_sequence = val_input_sequence;
config.val_output_sequence = val_output_sequence;
config.test_input_sequence = test_input_sequence;
config.test_output_sequence = test_output_sequence;
config.wash_out = wash_out;
config.err_type = err_type;
config.queue_type = queue_type;
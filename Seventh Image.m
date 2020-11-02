%% Testing Machine Learning in Matlab
clc; clearvars; close all;
%% Build Thermo Data sets

% training_size = 1e4;                    %number of sequences
% validation_size = 50;
% min_length = 15;                        %nt
% max_length = 15;                        %nt
% training_list = cell(training_size,5);  %preallocate
% validation_list = cell(validation_size,5);  %preallocate
%%
% tic
% for i = 1:training_size
%     currseq = randSeqFloat(randi([min_length,max_length]));
%     training_list{i,1}=currseq;
%     [training_list{i,2}, training_list{i,3}, training_list{i,4}]...
%         =AllawiThermoFloat(currseq);
%     training_list{i,5}=meltCalcIDTFloat(currseq,...
%         [training_list{i,2}, training_list{i,3}, training_list{i,4}],...
%         1e-4,1e-3,1.25e-2);
% end
% toc
% for i = 1:validation_size
%     currseq = randSeqFloat(randi([min_length,max_length]));
%     validation_list{i,1}=currseq;
%     [validation_list{i,2}, validation_list{i,3}, validation_list{i,4}]...
%         =AllawiThermoFloat(currseq);
%     validation_list{i,5}=meltCalcIDTFloat(currseq,...
%         [validation_list{i,2}, validation_list{i,3}, validation_list{i,4}],...
%         1e-4,1e-3,1.25e-2);
% end
% toc

%% Using double format
% seq_list_train = zeros(training_size,15);
% outputs = zeros(training_size,4);
% for i = 1:training_size
%     currseq = randSeqFloat(randi([min_length,max_length]));
%     seq_list_train(i,:)=currseq;
%     [outputs(i,1), outputs(i,2), outputs(i,3)]...
%         =AllawiThermoFloat(currseq);
%     outputs(i,4)=meltCalcIDTFloat(currseq,...
%         [outputs(i,1), outputs(i,2), outputs(i,3)],...
%         1e-4,1e-3,1.25e-2);
% end
% %results
% disp('NN prediction')
% disp(NN_Tm_test1(seq_list_train(1,:)))
% disp('Real vals')
% disp(outputs(1,:))
%% One time test of AllawiChar vs. AllawiFloat
% % sequence_list_char = cell(training_size,4);  %preallocate
% % for i = 1:training_size
% %     currseq = oligoCharConvert(sequence_list{i,1});
% %     sequence_list_char{i,1}=currseq;
% %     [sequence_list_char{i,2},sequence_list_char{i,3},sequence_list_char{i,4}]...
% %         =AllawiThermo(currseq);
% % end
%% Use zeros to indicate empty (terminal) nucleotides to give wider sequence arrays
training_size = 1e4;                    %number of sequences
min_length = 3;                        %nt
max_length = 30;                        %nt
seq_list_train = zeros(training_size,max_length);
outputs = zeros(training_size,4);

for i = 1:training_size
    randL = randi([min_length,max_length]);    %different each time
    currseq = randSeqFloat(randL);             %generate numeric sequence
    seq_list_train(i,1:randL)=currseq;         %slice is length of seq in sparse matrix
    [outputs(i,1), outputs(i,2), outputs(i,3)]...
        =AllawiThermoFloat(currseq);
    outputs(i,4)=meltCalcIDTFloat(currseq,...
        [outputs(i,1), outputs(i,2), outputs(i,3)],...
        1e-4,1e-3,1.25e-2);
end
%% test thermo prediction of NN2 on new sequence
randL = randi([min_length,max_length]);
currseq = randSeqFloat(randL);
slice = zeros(1,max_length);
seqwterminal0 = [currseq, zeros(1,max_length-randL)];
[outputsTest(1,1), outputsTest(1,2), outputsTest(1,3)]...
        =AllawiThermoFloat(currseq);
outputsTest(1,4)=meltCalcIDTFloat(currseq,...
        [outputsTest(1,1), outputsTest(1,2), outputsTest(1,3)],...
        1e-4,1e-3,1.25e-2);
disp(oligoCharConvert(currseq)')
disp('NN prediction')
disp(NN_Tm_test2(seqwterminal0))        %NN obtained from the NNfit app
disp('Real vals')
disp(outputsTest)
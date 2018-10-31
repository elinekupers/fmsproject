%% CreateFullFieldstim.m

% Change old stimulus presentation matrix to new one that only presents full field stim

for ii = 1:10
    
    load(fullfile(fmsRootPath, 'stimulus', '2014', sprintf('onOffLeftRight_params%d.mat',ii)));
    
    stimulus.seq(stimulus.seq == 5) = 1;
    stimulus.seq(stimulus.seq == 6) = 2;
    stimulus.seq(stimulus.seq == 7) = 1;
    stimulus.seq(stimulus.seq == 8) = 2;
    
    stimulus.trigSeq(stimulus.trigSeq == 5) = 1;
    stimulus.trigSeq(stimulus.trigSeq == 6) = 2;
    stimulus.trigSeq(stimulus.trigSeq == 7) = 1;
    stimulus.trigSeq(stimulus.trigSeq == 8) = 2;
    
    diodeSeq = zeros(size(stimulus.trigSeq));
    diodeSeq(stimulus.trigSeq == 1) = 1;
    diodeSeq(stimulus.trigSeq == 2) = 0;
    
    stimulus.diodeSeq = diodeSeq;
    
    save(fullfile(fmsRootPath, 'stimulus','2018', sprintf('onOffOnly_params%d.mat', ii)), 'stimulus');
    clear stimulus;
end
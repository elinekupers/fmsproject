function [subjectIDs, whichSession] = getSubjectIDs()

subjectIDs      = {'wlsubj002', ... % S1 - Full, Left, Right stim experiment
                   'wlsubj004', ... % S2 - Full, Left, Right stim experiment
                   'wlsubj005', ... % S3 - Full, Left, Right stim experiment
                   'wlsubj006', ... % S4 - Full, Left, Right stim experiment
                   'wlsubj010', ... % S5 - Full, Left, Right stim experiment
                   'wlsubj011', ... % S6 - Full, Left, Right stim experiment
                   'wlsubj048', ... % S7 - Full  stim only experiment
                   'wlsubj046', ... % S8 - Full  stim only experiment
                   'wlsubj039', ... % S9 - Full  stim only experiment
                   'wlsubj059', ... % S10 - Full  stim only experiment
                   'wlsubj067', ... % S11 - Full  stim only experiment
                   'wlsubj070'};    % S12 - Full  stim only experiment

% Go from subject to session nr
whichSession = [2, ... 'wlsubj002'
                7, ... 'wlsubj004'
                8, ... 'wlsubj005'
                1, ... 'wlsubj006'
                6, ... 'wlsubj010'
                5, ... 'wlsubj011'
                9, ... 'wlsubj048'
                10, ... 'wlsubj046'  % Full field Only
                11, ... 'wlsubj039'  % Full field Only
                12, ... 'wlsubj059'  % Full field Only
                13, ... 'wlsubj067'  % Full field Only
                14]; ...'wlsubj070'  % Full field Only
                
end
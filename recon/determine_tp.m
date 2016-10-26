    if exist('save_dir') & isdir(save_dir)
        for tp_idx = 1:ntps % Determine which tp to do next (check them all, every time).
            tp = tps_to_do(tp_idx);
            fn = dir([save_dir '/recon_tp' num2str(tp, '%.3d') '*.mat']);
            if length(fn) == 0 % This tp has not been done.
                save([save_dir '/recon_tp' num2str(tp, '%.3d') '.mat'], 'empty') % Create empty .mat file to prevent other instances of this script from computing the same tp.
                break; % Break out of for-loop and do this tp.
            end
            if tp_idx == ntps % All tps have been done.
                done_recon = 1;
            end
        end
        if done_recon
            break; % Break out of outer "while true" loop.
        end
    else % Don't use files in save_dir to determine which tp to do next.
        tp_idx = tp_idx + 1;
        if tp_idx > ntps
            break; % Break out of outer "while true" loop.
        end
        tp = tps_to_do(tp_idx);
    end

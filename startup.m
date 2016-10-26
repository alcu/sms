%%% Sets up irt toolbox.
orig_dir = pwd;
cd irt
setup
cd(orig_dir)
clear

%%% Adds some paths.
orig_dir = pwd;
addpath([orig_dir filesep 'tools'])

clear

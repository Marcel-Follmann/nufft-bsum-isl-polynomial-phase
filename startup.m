close all; clearvars; clc; 

addpath(genpath(pwd));
git_repo_check();
clean_asv_data();

function git_repo_check()
    % Get the status of the repository
    [~, statusMsg] = system('git status --porcelain');
    [~, summaryMsg] = system('git status -sb');
    
    % Check for untracked or uncommitted files
    fprintf("Welcome to this repository!\nPlease make sure to read the README before working with the code.\n\n");

    if isempty(strtrim(statusMsg))
        disp('✔️  No uncommitted or untracked changes.');
    else
        disp('❗️ There are untracked or uncommitted files. You should commit them!');
        return;
    end
    
    % Check if repo is ahead or behind remote
    if contains(summaryMsg, '[ahead')
        disp('⬆️  You have local commits that are not pushed yet. Please run "git push".');
    elseif contains(summaryMsg, '[behind')
        disp('⬇️  The remote has new commits! Please run "git pull".');
    else
        if isempty(strtrim(statusMsg))
            disp('✅ The repository is up to date!');
        end
    end
end

function clean_asv_data()
    files = dir(fullfile(pwd, '**', '*.asv'));
    for k = 1:length(files)
        delete(fullfile(files(k).folder, files(k).name));
    end
end
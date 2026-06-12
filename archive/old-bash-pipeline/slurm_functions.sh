alias list_all_jobs='sacct -u $(whoami) --format=JobID,JobName%100,User,State,Elapsed,AllocCPUS,NodeList,QOS,Start,End,MaxRSS --noheader | grep -v -E "batch|extern"'
alias list_current_jobs='squeue -u $(whoami) -o "%.18i %.50j %.10u %.8T %.10M %.6D %.10Q %.10q %.10m %R"'

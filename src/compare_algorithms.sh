./analyze_runs.sh
time matlab -nojvm -nodisplay -nosplash -singleCompThread -r "run $PWD/compare_algorithms.m" || { echo 'compare_algorithms failed' ; exit 1; }

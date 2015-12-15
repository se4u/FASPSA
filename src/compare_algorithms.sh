export dim_p=${1-120}
echo 'compare_algorithms.sh dim_p=' $dim_p
time matlab -nodisplay -nosplash -singleCompThread -r "dim_p = ${dim_p}; run $PWD/compare_algorithms.m" || { echo 'compare_algorithms failed' ; exit 1; }
./analyze_runs.sh $dim_p

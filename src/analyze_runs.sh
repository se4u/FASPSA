export dim_p=${1-15}
echo 'analyze_runs.sh dim_p=' $dim_p
matlab -singleCompThread -nosplash -r "dim_p=${dim_p}; run $PWD/analyze_runs.m"

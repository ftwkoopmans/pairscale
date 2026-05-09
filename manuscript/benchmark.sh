#!/bin/zsh

run_and_log() {
    # print the arguments separated by tabs (without a trailing newline)
    echo -n "$1\t$2\t$3\t" >> ~/temp/log_ram.txt

    # run the R script with time, pipe both stdout and stderr (|&) to grep and sed, and append to log
    /usr/bin/time -l Rscript benchmark.R "$1" "$2" "$3" |& grep peak | sed -E "s/^[^0-9]*([0-9]+)[^0-9]*peak.*/\1/" >> ~/temp/log_ram.txt
}

# iterate over all parameter combinations, repeating everything 3 times
for iter in 1 2 3; do
    for param_ncol in 100 500 1000 2000 4000; do
        for param_nrow in 10000 100000; do
            for param_algo in none msdap_mwmb vsn limma_loess pairscale_mean pairscale_median pairscale_trimmedmean pairscale_mode_fastest pairscale_mode_fast pairscale_mode_default; do

                # msdap_mwmb is slow; only apply it to small datasets
                if [[ "$param_algo" != "msdap_mwmb" ]] || (( param_ncol <= 2000 && param_nrow <= 10000 )); then
                    run_and_log "$param_algo" "$param_nrow" "$param_ncol"
                fi

            done
        done
    done
done

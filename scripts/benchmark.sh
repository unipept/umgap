#!/bin/bash


outdir="benchmark_results_$(date -Is)"
mkdir "$outdir"
pushd "$outdir"

u1=umgap-0.1.0  # git checkout 0.1.0
uio=umgap-0.1.1 # git checkout 0.1.1
u3=umgap-0.3.0  # git checkout 0.2.0

nc=openbsd-nc

#felixdir="/data/felix/metabenchmark"
#taxons="$felixdir/unipept-db/2018-06-taxons.tsv"
#tryptic_index="$felixdir/unipept-db/2018-06-tryptic.fst"
#ninemer_index_disk="$felixdir/unipept-db/2018-06-9mer.fst"
#ninemer_index="/data/tempfs/2018-06-9mer.fst"

felixdir="/data/felix/metabenchmark"
taxons="/home/rien/Unipept/data/2018-06-taxons.tsv"
tryptic_index="/home/rien/Unipept/data/2018-06-tryptic.fst"
ninemer_index_disk="/home/rien/Unipept/data/2018-06-tryptic.fst"
ninemer_index_tempfs="/tmp/2018-06-tryptic.fst"

#data_pe1="/home/rien/Unipept/data/mini_1.fq"
#data_pe2="/home/rien/Unipept/data/mini_2.fq"
data_pe1="/home/rien/Unipept/data/small_1.fq"
data_pe2="/home/rien/Unipept/data/small_2.fq"

full_analysis_log="full_analysis_log.tsv"

ncpus="$(cat /proc/cpuinfo  | grep processor | wc -l)"
maxcpu="$(($ncpus-1))"

echo "[$(date)] Benchmarking with $ncpus cpu's"

# Full pipeline analysis

umgap1_tempfs() {
    $u1 fastq2fasta "$data_pe1" "$data_pe2"         | # intercalate paired ends
      sed '/^>/s_/.*__'                             | # remove paired end indicator
      $u1 translate -a                              | # translate 6 frames
      $u1 prot2kmer2lca -o -k9 "$ninemer_index_tempfs"  | # map to taxa
      $u1 seedextend -g4 -s4                        | # filter with seed extend
      $u1 uniq                                      | # join ends/frames
      $u1 taxa2agg -l2 -a lca\* "$taxons"           | # aggregate
      grep -v '^>'                                  | # remove headers
      $u1 report "$taxons"                            # create report
}

umgap3_tempfs() {
    $u3 fastq2fasta "$data_pe1" "$data_pe2"         | # intercalate paired ends
      sed '/^>/s_/.*__'                             | # remove paired end indicator
      $u3 translate -a                              | # translate 6 frames
      $u3 prot2kmer2lca -o -k9 "$ninemer_index_tempfs"  | # map to taxa
      $u3 seedextend -g4 -s4                        | # filter with seed extend
      $u3 uniq                                      | # join ends/frames
      $u3 taxa2agg -l2 -a lca\* "$taxons"           | # aggregate
      grep -v '^>'                                  | # remove headers
      $u3 report "$taxons"                            # create report
}

umgap3_inmemory() {
    $u3 fastq2fasta "$data_pe1" "$data_pe2"         | # intercalate paired ends
      sed '/^>/s_/.*__'                             | # remove paired end indicator
      $u3 translate -a                              | # translate 6 frames
      $u3 prot2kmer2lca -m -o -k9 "$ninemer_index_disk"  | # map to taxa
      $u3 seedextend -g4 -s4                        | # filter with seed extend
      $u3 uniq                                      | # join ends/frames
      $u3 taxa2agg -l2 -a lca\* "$taxons"           | # aggregate
      grep -v '^>'                                  | # remove headers
      $u3 report "$taxons"                            # create report
}

umgap3_socket() {
    $u3 fastq2fasta "$data_pe1" "$data_pe2"         | # intercalate paired ends
      sed '/^>/s_/.*__'                             | # remove paired end indicator
      $u3 translate -a                              | # translate 6 frames
      $nc -NU umgap-socket                          | # map to taxa
      $u3 seedextend -g4 -s4                        | # filter with seed extend
      $u3 uniq                                      | # join ends/frames
      $u3 taxa2agg -l2 -a lca\* "$taxons"           | # aggregate
      grep -v '^>'                                  | # remove headers
      $u3 report "$taxons"                            # create report
}

# Pipeline subcommands

fastq2fasta() {
    $umgap fastq2fasta "$data_pe1" "$data_pe2" | sed '/^>/s_/.*__'
}

translate(){
    $umgap translate -a
}

prot2kmer2lca_tempfs(){
    $umgap prot2kmer2lca -o -k9 "$ninemer_index_tempfs"
}

prot2kmer2lca_inmemory(){
    $umgap prot2kmer2lca -m -o -k9 "$ninemer_index_disk"
}

prot2kmer2lca_socket(){
      $nc -NU umgap-socket
}

seedextend(){
    $umgap seedextend -g4 -s4
}

uniq(){
    $umgap uniq
}

taxa2agg(){
    $umgap taxa2agg -l2 -a lca\* "$taxons"
}

report(){
    grep -v '^>' | $umgap report "$taxons"
}

# Helper functions

sort_fasta(){
    sed '/^\s*$/d' "$1"             | # Remove empty lines
        sed -z "s/\n\([^>]\)/~\1/g" | # Put sequences on the same line as their header
        sort                        | # Sort
        sed "s/~/\n/"                 # Restore sequences
}

set_cpus(){
    enabled_cpus="$1"
    # current PID and all background PID's
    for pid in $$ $(jobs -lp | tr '\n' ' '); do
        taskset -apc "$enabled_cpus" "$pid"
    done
}

time_part(){
    cmd="$1"
    input="$2"
    output="$3"
    echo "[$(date)] $cmd starting"
    TIMEFORMAT="$1	%R"
    time cat "$input" | $cmd > "$output"
    echo "[$(date)] $cmd done"
}

time_tool(){
    echo "[$(date)] $1 starting"
    TIMEFORMAT="%R"
    time $@
    echo "[$(date)] $1 done"
}

setup_tempfs() {
    echo "[$(date)] Copying to tempfs"
    #mount /data/tempfs
    cp "$ninemer_index_disk" "$ninemer_index_tempfs"
    echo "[$(date)] Copied to tempfs"
}

breakdown_tempfs() {
    echo "[$(date)] Removing from tempfs"
    rm "$ninemer_index_tempfs"
    #umount /data/tempfs
    echo "[$(date)] Removed from tempfs"
}


compare_results(){
    for i in $(seq 1 6); do
        if ! cmp <(sort_fasta "$1_part_$i") <(sort_fasta "$2_part_$i"); then
            echo "[$(date)] WARNING: Intermediate results (part_$i) differ between $1 and $2"
            cat "$1_part_$i"
            echo "~~~"
            cat "$2_part_$i"
        fi
    done

    if ! cmp "$1_report" "$2_report"; then
        echo "[$(date)] WARNING: Final reports differ between $1 and $2"
            cat "$1_report"
            echo "~~~"
            cat "$2_report"
    fi
}

# Benchmark functions

benchmark_each_part(){
    umgap="$1"
    base="$2"
    logfile="${base}_time.tsv"
    echo "[$(date)] Starting benchmark '$base'"
    echo -e "command\ttime(s)" >"$logfile"
    time_part fastq2fasta           /dev/null          "${base}_part_1"   2>>"$logfile"
    time_part translate             "${base}_part_1"   "${base}_part_2"   2>>"$logfile"
    time_part prot2kmer2lca_tempfs  "${base}_part_2"   "${base}_part_3"   2>>"$logfile"
    time_part seedextend            "${base}_part_3"   "${base}_part_4"   2>>"$logfile"
    time_part uniq                  "${base}_part_4"   "${base}_part_5"   2>>"$logfile"
    time_part taxa2agg              "${base}_part_5"   "${base}_part_6"   2>>"$logfile"
    time_part report                "${base}_part_6"   "${base}_report"   2>>"$logfile"
}

benchmark_strong_scaling() {
    tool="$1"
    in="$2"
    base="$3"
    logfile="${base}_strong_scaling.tsv"
    out="/dev/null"

    echo "[$(date)] Starting strong scaling benchmark '$base'"

    echo -e "cpus\ttime(s)" > "$logfile"

    cpus=()
    for i in $(seq 0 $maxcpu); do
        cpus+=($i)
        enabled_cpus=$(IFS=, ; echo "${cpus[*]}")
        echo "[$(date)] Strong scaling of $tool with cpus: $enabled_cpus"
        set_cpus "$enabled_cpus"
        echo -en "$i\t" >> "$logfile"
        cat $input | time_tool "$tool" "$in" >/dev/null 2>> "$logfile"
    done
    set_cpus "0-$maxcpu"
}

benchmark_weak_scaling() {
    tool="$1"
    input_full="$2"
    base="$3"

    echo "Starting weak scaling benchmark '$base'"

    logfile="${base}_weak_scaling.tsv"
    echo -e "ncpus($lines_per_cpu lines per cpu)\ttime(s)" > $logfile

    input_lines="$(wc -l $input_full | cut -d' ' -f1)"
    lines_per_cpu="$(($input_lines/$ncpus))"
    if [[ ! $(($input_lines % $ncpus)) -eq 0 ]]; then
        echo "$input_full is not divisabe in equal parts with $ncpus cpu's"
        exit 1
    fi

    echo -e "cpus\ttime(s)" > "$logfile"

    cpus=()
    for i in $(seq 0 $maxcpu); do
        cpus+=($i)
        enabled_cpus=$(IFS=, ; echo "${cpus[*]}")

        end=$(($lines_per_cpu*($i + 1)))
        sed -n "1,${end}p" $input_full > partial_input

        echo "[$(date)] Weak scaling of $tool with cpus $enabled_cpus and $end lines"
        set_cpus "$enabled_cpus"
        echo -en "$(($i+1))\t" >> "$logfile"
        cat $input | time_tool "$tool" "$in" >/dev/null 2>> "$logfile"
    done
    set_cpus "0-$maxcpu"
}

function benchmark_full_analysis(){
    tool="$1"
    echo "[$(date)] Performing a full analysis benchmark with '$1'"
    echo -en "$tool\t" >>$full_analysis_log
    time_tool $tool > /dev/null 2>>$full_analysis_log
}

function safe_exit(){
    breakdown_tempfs || true
    exit 1
}


trap safe_exit SIGINT

#
# Measure improvement per subcommand & check if the results are still the same
#

setup_tempfs

benchmark_each_part $u1  "umgap_v1_subcommands"
benchmark_each_part $uio "umgap_io_subcommands"
benchmark_each_part $u3  "umgap_v3_subcommands"

compare_results "umgap_v1_subcommands" "umgap_io_subcommands"
compare_results "umgap_v1_subcommands" "umgap_v3_subcommands"

dna_data="umgap_v1_subcommands_part_1"
protein_data="umgap_v1_subcommands_part_2"

# Measure parallel scaling of prot2kmer2lca

umgap=$u1
benchmark_strong_scaling "prot2kmer2lca_tempfs" $protein_data "umgap_v1"
benchmark_weak_scaling "prot2kmer2lca_tempfs" $protein_data "umgap_v1"

umgap=$u3
benchmark_strong_scaling "prot2kmer2lca_tempfs" $protein_data "umgap_v3"
benchmark_weak_scaling "prot2kmer2lca_tempfs" $protein_data "umgap_v3"

# Full benchmark

benchmark_full_analysis umgap1_tempfs
benchmark_full_analysis umgap3_tempfs

breakdown_tempfs

benchmark_full_analysis umgap3_inmemory

echo "[$(date)] Starting umgap daemon ..."

$u3 prot2kmer2lca -m -o -k9 --socket umgap-socket "$ninemer_index_disk" &

echo "[$(date)] Waiting untill socket is ready..."
while [ ! -S umgap-socket ]; do sleep 1; done
echo "[$(date)] Socket ready, benchmarking with socket..."

benchmark_full_analysis umgap3_socket

benchmark_strong_scaling prot2kmer2lca_socket $protein_data "umgap_v3_socket"
benchmark_weak_scaling prot2kmer2lca_socket $protein_data "umgap_v3_socket"

echo "[$(date)] Killing umgap daemon"

kill %1
rm umgap-socket

popd


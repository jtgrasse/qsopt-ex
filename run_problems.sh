#!/usr/bin/env bash
# Solve every problem in Problems/ with build/esolver/esolver and write, for
# each problem, to Problems_Outputs/:
#   <problem>_bases.txt  pivot log from esolver -W (A, b, c, bases, updates,
#                        final exact solution)
#   <problem>.log        everything esolver printed to the terminal
# One solve runs per L3 cache at a time, pinned to that cache's CPUs, so
# concurrent solves never share an L3.  Largest problems are started first.
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
esolver=$root/build/esolver/esolver
probdir=$root/Problems
outdir=$root/Problems_Outputs

if [[ ! -x $esolver ]]; then
	echo "error: $esolver not found, build the project in build/ first" >&2
	exit 1
fi
mkdir -p "$outdir"

# print one CPU list (e.g. 0-7) per distinct L3 cache
l3_groups() {
	local d
	for d in /sys/devices/system/cpu/cpu[0-9]*/cache/index*; do
		[[ -r $d/level && $(<"$d/level") == 3 ]] && cat "$d/shared_cpu_list"
	done | sort -u -n
}

mapfile -t groups < <(l3_groups)
if (( ${#groups[@]} == 0 )); then
	echo "warning: no L3 cache information found, running one solve unpinned" >&2
	groups=("")
fi

# problem name as esolver derives it: basename without .gz/.bz2 and extension
problem_name() {
	local n
	n=$(basename "$1")
	n=${n%.gz}
	n=${n%.bz2}
	echo "${n%.*}"
}

mapfile -t probs < <(find "$probdir" -maxdepth 1 -type f -printf '%s\t%p\n' |
	sort -rn | cut -f2-)
if (( ${#probs[@]} == 0 )); then
	echo "no problems found in $probdir" >&2
	exit 1
fi

# memory limit per solve (esolver -m, bytes): half of the total memory, split
# evenly between the concurrent workers
memtotal_kb=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)
memlimit=$((memtotal_kb * 1024 / 2 / ${#groups[@]}))

echo "${#probs[@]} problem(s), ${#groups[@]} L3 cache(s): ${groups[*]:-unpinned}"
echo "memory limit per solve: $((memlimit / 1024 / 1024)) MiB"

# shared work queue: index of the next problem, guarded by flock
queue=$(mktemp)
echo 0 > "$queue"
trap 'rm -f "$queue" "$queue.lock"' EXIT
trap 'echo "interrupted, stopping solves" >&2; kill 0' INT TERM

next_index() {
	local i
	{
		flock 9
		i=$(<"$queue")
		echo $((i + 1)) > "$queue"
	} 9>"$queue.lock"
	echo "$i"
}

# solve problems from the queue on the CPUs of one L3 cache
worker() {
	local cpus=$1 i f name rc start
	local pin=()
	[[ -n $cpus ]] && pin=(taskset -c "$cpus")
	while :; do
		i=$(next_index)
		(( i < ${#probs[@]} )) || break
		f=${probs[i]}
		name=$(problem_name "$f")
		echo "[cpus ${cpus:-all}] start  $name"
		start=$SECONDS
		rc=0
		(cd "$outdir" && "${pin[@]}" "$esolver" -m "$memlimit" -W "${name}_bases.txt" "$f") \
			> "$outdir/$name.log" 2>&1 || rc=$?
		echo "[cpus ${cpus:-all}] finish $name: exit $rc, $((SECONDS - start)) s"
	done
}

for g in "${groups[@]}"; do
	worker "$g" &
done
wait
echo "all done, outputs in $outdir"

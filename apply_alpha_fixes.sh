#!/usr/bin/env bash
# Fold the validation run's side-config workarounds into the pipeline itself,
# so the validation exercises the released code rather than local overrides.
#
# Three changes:
#   1. conf/modules.config          72 h -> 168 h on the four long-running tasks
#   2. modules/nf-core/iqtree       -nt AUTO -ntmax N  ->  -nt N
#   3. modules/nf-core/gubbins      export TMPDIR=/tmp before run_gubbins.py
#
# Run from the repository root. Every edit asserts on an exact match, so the
# script fails loudly rather than half-applying.

set -euo pipefail

if [ ! -f nextflow.config ] || [ ! -d modules/nf-core ]; then
    echo "ERROR: run this from the CholeraSeq2 repository root" >&2
    exit 1
fi

python3 - <<'PYEOF'
import sys

def edit(path, old, new, expect=1):
    with open(path) as fh:
        text = fh.read()
    n = text.count(old)
    if n != expect:
        sys.exit(f"ERROR: {path}: expected {expect} occurrence(s) of the target "
                 f"text, found {n}. Nothing written.")
    with open(path, 'w') as fh:
        fh.write(text.replace(old, new))
    print(f"  {path}: {n} edit(s)")


# ---------------------------------------------------------------- 1. wall time
# Gubbins on the 4,266-sequence cluster.1 spent 9 h 34 m on RAxML for iteration
# 1 alone at 6 threads, and 40 h at 5. With up to five iterations plus ancestral
# reconstruction, and nothing cached until the task completes, 72 h is too tight
# to be safe. One week on each of the four long-running tasks.
print("conf/modules.config")
edit('conf/modules.config', 'time = 72.h', 'time = 168.h', expect=4)
edit('conf/modules.config',
     'Allow 72 h (still capped by params.max_time).',
     'Allow 168 h (still capped by params.max_time).',
     expect=2)
edit('conf/modules.config',
     '// Give it 72 h by default;',
     '// Give it one week by default;',
     expect=1)


# -------------------------------------------------------------- 2. IQ-TREE -nt
# Under `-nt AUTO`, IQ-TREE sizes its thread pool to the cores it detects on the
# host before running the benchmark that would apply -ntmax, so the task is not
# bound by its cpus allocation. On this 32-core host a task given 6 cpus was
# measured at 1963% CPU across 32 threads, and spent over 46 h in a stage that
# completes in 52 s with an explicit -nt. An explicit -nt makes -ntmax redundant.
print("modules/nf-core/iqtree/main.nf")
edit('modules/nf-core/iqtree/main.nf',
     "        -nt AUTO \\\\\n        -ntmax $task.cpus \\\\\n",
     "        -nt $task.cpus \\\\\n")
edit('modules/nf-core/iqtree/main.nf',
     "    def memory                      = task.memory.toString().replaceAll(' ', '')\n",
     "    def memory                      = task.memory.toString().replaceAll(' ', '')\n"
     "    // The thread count is given explicitly rather than with `-nt AUTO`.\n"
     "    // Under AUTO, IQ-TREE sizes its thread pool to the cores it detects on\n"
     "    // the machine before the benchmark that would apply -ntmax, so the task\n"
     "    // is not bound by its cpus allocation: on a 32-core host a task\n"
     "    // allocated 6 cpus was measured at 1963% CPU across 32 threads, and\n"
     "    // spent over 46 h in a stage that takes 52 s with an explicit -nt.\n"
     "    // An explicit -nt also makes -ntmax redundant.\n")


# -------------------------------------------------------------- 3. Gubbins tmp
# pyjar runs under Python multiprocessing, whose Manager creates its temp
# directory under TMPDIR. Nextflow points TMPDIR at the task work directory;
# where that is on NFS, a file deleted while still open is silly-renamed to
# .nfsXXXXXX and cannot be unlinked, so the Manager's shutdown cleanup raises
# OSError: [Errno 16] Device or resource busy. The Manager then exits while its
# workers are still running and the task hangs indefinitely rather than failing,
# so no exit status is ever written and only the wall-clock limit ends it.
print("modules/nf-core/gubbins/main.nf")
edit('modules/nf-core/gubbins/main.nf',
     '    export PYTHONPATH="\\$PWD/pyfix\\${PYTHONPATH:+:\\$PYTHONPATH}"\n',
     '    export PYTHONPATH="\\$PWD/pyfix\\${PYTHONPATH:+:\\$PYTHONPATH}"\n'
     '\n'
     "    # pyjar's ancestral reconstruction runs under Python multiprocessing, whose\n"
     '    # Manager process creates its temp directory under TMPDIR. Nextflow points\n'
     '    # TMPDIR at the task work directory, and where that is a network filesystem a\n'
     '    # file deleted while still open is silly-renamed to .nfsXXXXXX and cannot be\n'
     '    # unlinked, so the Manager\'s shutdown cleanup raises\n'
     "    #   OSError: [Errno 16] Device or resource busy: '.nfsXXXXXX'\n"
     '    # The Manager then exits while its workers are still running and the task\n'
     '    # hangs indefinitely instead of failing: no exit status is written, one core\n'
     '    # idles, and only the wall-clock limit ever ends it. Keeping the Manager\'s\n'
     '    # temp directory on container-local storage avoids the rename entirely.\n'
     '    export TMPDIR=/tmp\n')

print("\nAll edits applied.")
PYEOF

echo
echo "Review before committing:"
echo "  git diff --stat"
echo "  git diff"

#! /bin/shA

source $HOME/local/ase-svn/env.sh

ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

PYTHONPATH="$PYTHONPATH:$ROOT"

# main.py

# usage pattern
usage = '''

Usage:
  docopt_example.py command --option <argument>
  docopt_example.py <argument> <repeating-argument>
  docopt_example.py --version

Options:
  -h, --help     Display help
  -o, --option   Display options
  -l, --all      List all
  -q, --quit     exit
  --version      Version 3.6.1

'''

# Initialization
from docopt import docopt

args = docopt(usage)
print(args)
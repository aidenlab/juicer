#!/usr/bin/env python

# Generate site positions in genome from given restriction enzyme
# Juicer 1.5

from __future__ import print_function

import sys
import re

def usage():
  print('Usage: {} <restriction enzyme> <genome> [location]'.format(sys.argv[0]), file=sys.stderr)
  sys.exit(1)

# ------------------------------------------------------------------------------

def process_args(args):

  # Genome to filename mappings
  #
  # You may hardcode filenames belonging to frequently used genomes by inserting
  # elements into this dictionary.

  filenames = {
    'hg19': '/seq/references/Homo_sapiens_assembly19.fasta',
    'mm9' : '/seq/references/Mus_musculus_assembly9.fasta',
    'mm10': '/seq/references/Mus_musculus_assembly10.fasta',
    'hg18': '/seq/references/Homo_sapiens_assembly18.fasta',
  }

  # Enzyme to search pattern mappings
  #
  # You may provide your own patterns by inserting elements into this dictionary.
  # Multiple patterns are supported in the form of lists. If you enumerate more
  # than one patterns for an enzyme, then a match will be reported when at least
  # one of them can be found at a given position.

  patterns = {
    'HindIII'  : 'AAGCTT',
    'DpnII'    : 'GATC',
    'MboI'     : 'GATC',
    'Sau3AI'   : 'GATC',
    'MultiTest': [ 'AAGCTT', 'GATC' ],
  }

  if len(args) != 3 and len(args) != 4:
    usage()

  enzyme     = args[1]
  genome     = args[2]
  inputfile  = ''
  outputfile = ''
  pattern    = ''

  if len(args) == 4:
    inputfile = args[3]
  elif genome in filenames:
    inputfile = filenames[genome]
  else:
    print('<genome> not found and [location] not defined', file=sys.stderr)
    usage()

  if enzyme in patterns:
    if isinstance(patterns[enzyme], str):
      pattern = patterns[enzyme].lower()
    else:
      pattern = [ p.lower() for p in set(patterns[enzyme]) ]
  else:
    print('<restriction enzyme> must be one of {}'.format(list(patterns.keys())), file=sys.stderr)
    usage()

  outputfile = genome + '_' + enzyme + '.txt'

  return {
    'enzyme'     : enzyme,
    'genome'     : genome,
    'pattern'    : pattern,
    'inputfile'  : inputfile,
    'outputfile' : outputfile,
  }

# ------------------------------------------------------------------------------

def process_single_pattern(params):

  f = open(params['inputfile' ], 'r')
  g = open(params['outputfile'], 'w')

  pattern    = params['pattern']
  patternlen = len(pattern)

  segment = ''
  counter = 0
  endl    = ''

  for line in f:

    line = line.strip()

    if line.startswith('>'):

      # This is the start of a new sequence.

      if counter > 0:
        g.write(' ' + str(counter)) # Close the previous sequence here.

      g.write(endl+line[1:]) # Start the new sequence here.

      segment = ''
      counter = 0
      endl    = '\n'

      continue

    # Process next line of the sequence.

    line = line.lower()

    for symbol in line:

      counter += 1
      segment += symbol

      while len(segment) > patternlen:
        segment = segment[1:]
      if segment == pattern:
        g.write(' ' + str(counter - patternlen + 1)) # There is a match here.

  # Close the last sequence.

  if counter > 0:
    g.write(' ' + str(counter))

  g.write('\n') # End the output file with a newline.

  # Close files.

  g.close()
  f.close()

# ------------------------------------------------------------------------------

def process_multi_pattern(params):

  f = open(params['inputfile' ], 'r')
  g = open(params['outputfile'], 'w')

  patterns = params['pattern']
  maxsize  = 0                 # Maximum segment size (size of the longest pattern).
  minsize  = len(patterns[0])  # Minimum segment size (size of the shortest pattern).

  for p in patterns:
    if len(p) > maxsize:
      maxsize = len(p)
    if len(p) < minsize:
      minsize = len(p)

  segment = ''
  counter = 0
  endl    = ''

  for line in f:

    line = line.strip()

    if line.startswith('>'):

      # This is the beginning of a new sequence, but before starting it we must
      # finish processing of the remaining segment of the previous sequence.

      while len(segment) > minsize:
        segment = segment[1:]
        match   = False
        for p in patterns:
          if segment.startswith(p):
            match = True
            break
        if match:
          g.write(' ' + str(counter - len(segment) + 1))

      if counter > 0:
        g.write(' ' + str(counter)) # Close the previous sequence here.

      g.write(endl+line[1:]) # Start the new sequence here.

      segment = ''
      counter = 0
      endl    = '\n'

      continue

    # Process next line of the sequence.

    line = line.lower()

    for symbol in line:

      counter += 1
      segment += symbol

      while len(segment) > maxsize:
        segment = segment[1:]

      # Do pattern matching only if segment size equals maxsize.

      if len(segment) == maxsize:
        match = False
        for p in patterns:
          if segment.startswith(p):
            match = True
            break
        if match:
          g.write(' ' + str(counter - maxsize + 1)) # maxsize == len(segment)

  # Finish the last sequence.

  while len(segment) > minsize:
    segment = segment[1:]
    match   = False
    for p in patterns:
      if segment.startswith(p):
        match = True
        break
    if match:
      g.write(' ' + str(counter - len(segment) + 1))

  if counter > 0:
    g.write(' ' + str(counter))

  g.write('\n') # End the output file with a newline.

  # Close files.

  g.close()
  f.close()

# ------------------------------------------------------------------------------

params = process_args(sys.argv)

if isinstance(params['pattern'], str):
  process_single_pattern(params)
else:
  process_multi_pattern(params)

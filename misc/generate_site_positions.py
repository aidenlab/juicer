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

def pattern2regexp(pattern):

  wildcards = {
    'N': '[ACGT]',
    'M': '[AC]',
    'R': '[AG]',
    'W': '[AT]',
    'Y': '[CT]',
    'S': '[CG]',
    'K': '[GT]',
    'H': '[ACT]',
    'B': '[CGT]',
    'V': '[ACG]',
    'D': '[AGT]',
  }

  if isinstance(pattern, str):
    plen    = len(pattern)
    pattern = pattern.upper()
    for p, r in wildcards.items():
      pattern = re.sub(p, r, pattern)
    # Create a tuple containing the regular expression and the original pattern length.
    return (re.compile(pattern.lower()), plen)

  return [ pattern2regexp(p) for p in pattern ]

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
  #
  # Wildcards:
  #   N: A, C, G or T
  #   M: A or C
  #   R: A or G
  #   W: A or T
  #   Y: C or T
  #   S: C or G
  #   K: G or T
  #   H: A, C or T
  #   B: C, G or T
  #   V: A, C or G
  #   D: A, G or T

  patterns = {
    'HindIII'     : 'AAGCTT',
    'DpnII'       : 'GATC',
    'MboI'        : 'GATC',
    'Sau3AI'      : 'GATC',
    'Arima'       : [ 'GATC', 'GANTC' ],
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
    pattern = pattern2regexp(patterns[enzyme])
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

  pattern    = params['pattern'][0] # A regular expression.
  patternlen = params['pattern'][1] # An integer number.

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

      if re.match(pattern, segment):
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

  patterns = [ p[0] for p in params['pattern'] ]
  maxsize  = 0                       # Maximum segment size (size of the longest pattern).
  minsize  = params['pattern'][0][1] # Minimum segment size (size of the shortest pattern).

  for p in params['pattern']:
    if p[1] > maxsize:
      maxsize = p[1]
    if p[1] < minsize:
      minsize = p[1]

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
        for p in patterns:
          if re.match(p, segment):
            g.write(' ' + str(counter - len(segment) + 1))
            break

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
        for p in patterns:
          if re.match(p, segment):
            g.write(' ' + str(counter - maxsize + 1)) # maxsize == len(segment)
            break

  # Finish the last sequence.

  while len(segment) > minsize:
    segment = segment[1:]
    for p in patterns:
      if re.match(p, segment):
        g.write(' ' + str(counter - len(segment) + 1))
        break

  if counter > 0:
    g.write(' ' + str(counter))

  g.write('\n') # End the output file with a newline.

  # Close files.

  g.close()
  f.close()

# ------------------------------------------------------------------------------

params = process_args(sys.argv)

if isinstance(params['pattern'], list):
  process_multi_pattern(params)
else:
  process_single_pattern(params)

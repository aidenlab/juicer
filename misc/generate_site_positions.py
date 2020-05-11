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
    pattern = patterns[enzyme]
    # Convert a simple string to a list.
    if not isinstance(pattern, list):
      pattern = [ pattern ]
    # Make patterns uppercase.
    pattern = [ p.upper() for p in pattern ]
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

def has_wildcard(pattern):

  # Input pattern can be a list or a string.

  wildcards = re.compile(r'[NMRWYSKHBVD]')

  if (isinstance(pattern, list)):
    for p in pattern:
      if re.search(wildcards, p):
        return True
  else:
    if re.search(wildcards, pattern):
      return True

  return False

# ------------------------------------------------------------------------------

def pattern2regexp(pattern):

  # Input pattern can be a list or a string.

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

  if isinstance(pattern, list):
    return [ pattern2regexp(p) for p in pattern ]

  pattern = pattern.upper()
  for p, r in wildcards.items():
    pattern = re.sub(p, r, pattern)

  return re.compile(pattern.upper())

# ------------------------------------------------------------------------------

def get_match_func(pattern):

  # Input pattern can be a list or a string.

  if not isinstance(pattern, list):
    pattern = [ pattern ]

  if has_wildcard(pattern):

    pattern = pattern2regexp(pattern)

    if len(pattern) == 1: # There is only a single pattern.

      pattern = pattern[0] # Use the only element from the list as a single regexp.

      def match_single_regexp(segment):
        if re.match(pattern, segment):
          return True
        return False

      return match_single_regexp

    else: # There are multiple patterns.

      def match_multi_regexp(segment):
        for p in pattern:
          if re.match(p, segment):
            return True
        return False

      return match_multi_regexp

  else: # No wildcard in any of the patterns.

    if len(pattern) == 1: # There is only a single pattern.

      pattern = pattern[0] # Use the only element from the list as a single string.

      def match_single_string(segment):
        if segment.startswith(pattern):
          return True
        return False

      return match_single_string

    else: # There are multiple patterns.

      def match_multi_string(segment):
        for p in pattern:
          if segment.startswith(p):
            return True
        return False

      return match_multi_string

# ------------------------------------------------------------------------------

def process_input(params):

  f = open(params['inputfile' ], 'r')
  g = open(params['outputfile'], 'w')

  minsize = min([ len(p) for p in params['pattern'] ])
  maxsize = max([ len(p) for p in params['pattern'] ])
  matches = get_match_func(params['pattern'])

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
        if matches(segment):
          g.write(' ' + str(counter - len(segment) + 1))

      if counter > 0:
        g.write(' ' + str(counter)) # Close the previous sequence here.

      firststr=re.split('\s+',line[1:])
      g.write(endl+firststr[0])

      segment = ''
      counter = 0
      endl    = '\n'

      continue

    # Process next line of the sequence.

    line = line.upper()

    for symbol in line:

      counter += 1
      segment += symbol

      while len(segment) > maxsize:
        segment = segment[1:]

      # Do pattern matching only if segment size equals maxsize.

      if len(segment) == maxsize:
        if matches(segment):
          g.write(' ' + str(counter - maxsize + 1)) # maxsize == len(segment)

  # Finish the last sequence.

  while len(segment) > minsize:
    segment = segment[1:]
    if matches(segment):
      g.write(' ' + str(counter - len(segment) + 1))

  if counter > 0:
    g.write(' ' + str(counter))

  g.write('\n') # End the output file with a newline.

  # Close files.

  g.close()
  f.close()

# ------------------------------------------------------------------------------

params = process_args(sys.argv)
process_input(params)

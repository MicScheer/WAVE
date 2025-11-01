#     Utility to write the current results of wave.adjust to wave_adjust.dat

import sys

def wave_adjust(argv):

  fin = "wave.out" # default file to read

  nargs = len(argv)
  if nargs == 2: fin = argv[1] # file to read

  ffin = open(fin,'r')
  lines = ffin.readlines()
  ffin.close()

  bymax = 9999.
  bzmax = 9999.

  for line in lines:

      chkey = "BYmax,BYmin:"
      i1 = line.find(chkey)
      if i1 > -1:
        words = line.split(':')
        vals = words[1].split()
        bymax = float(vals[0])
        bymin = float(vals[1])
        continue

      chkey = "BZmax,BZmin:"
      i1 = line.find(chkey)
      if i1 > -1:
        words = line.split(':')
        vals = words[1].split()
        bzmax = float(vals[0])
        bzmin = float(vals[1])
        continue

  #endfor line in lines

# This chi2 is minimized by WAVE

  chi2 = ((bymax-bzmax)/(bymax+bzmax))**2
  print("\nBYmax, BZmax, Chi2:",bymax,bzmax,chi2)

  fchi2 = open('wave_adjust.dat','w')
  fchi2.write(str(chi2)+'\n')
  fchi2.close()

#enddef wave_adjust(argv)

wave_adjust(sys.argv)

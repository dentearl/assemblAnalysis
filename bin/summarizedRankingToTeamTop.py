#!/usr/bin/env python
"""
"""
import sys
import signal # deal with broken pipes

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def main():
   prevLineEmpty = False
   for line in sys.stdin:
      line = line.strip()
      if line.startswith('#'):
         print line
      if line == '':
         prevLineEmpty = True
         continue
      if prevLineEmpty == True:
         print line
         prevLineEmpty = False

if __name__ == '__main__':
   main()

import argparse
import MSA


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional file")
parser.add_argument("-m", "--miranda", help="Miranda file")
parser.add_argument("-t", "--ts", help="TargetSpy file")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()


def main():
  pc = MSA.PositionalConsensus(args.file)
  pc.get_full_target_region(args.miranda, args.ts, args.output)

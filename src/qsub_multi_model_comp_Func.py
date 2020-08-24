from multi_model_comp_Func import Model_Run
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("RUN_ID")
parser.add_argument("VALnits")
parser.add_argument("VALpits")
parser.add_argument("VALburnin")
args = parser.parse_args()

print('RUN_ID, VALnits, VALpits, VALburnin')
print(args.RUN_ID, args.VALnits, args.VALpits, args.VALburnin, sep=',')
#
Model_Run(str(args.RUN_ID), int(args.VALnits), int(args.VALpits), int(args.VALburnin))

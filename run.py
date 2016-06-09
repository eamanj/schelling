#import matplotlib.pyplot as plt
import argparse
import schelling

parser = argparse.ArgumentParser(description="")
parser.add_argument("-s" , "--size", dest="size", default=30, type=int,
                    help="size of the grid.")
parser.add_argument("-e" , "--empty_ratio", dest="empty_ratio", default=0.3, type=float,
                    help="Empty ratio of the grid")
parser.add_argument("-d", "--distribution", dest="distribution", default="Powerlaw",
                    choices=["Powerlaw", "Gaussian"], 
                    help="The distribution of initial similarit thresholds")
parser.add_argument("-a", "--alpha", dest="alpha", default=0.4, type=float,
                    help="alpha of power law")
parser.add_argument("-m", "--mean", dest="mean", default=0.4, type=float,
                    help="The mean of guassian distribution")
parser.add_argument("-st", "--std", dest="std", default=0.1, type=float,
                    help="The std of guassian distribution")
parser.add_argument("-i", "--max_iterations", dest="max_iters", default=600, type=int,
                    help="max number of iterations")
parser.add_argument("-r", "--num_races", dest="num_races", default=2,
                    help="number of races")
parser.add_argument("-u", "--update_rate", dest="update_rate", default=None, type=float,
                    help="Determines the update rate of tolerances. The default value is "
                    "None which indicates no updating (fixed tolerances). Should be a "
                    "float with between 0 and 1.") 
args = parser.parse_args()

def main():
  model = schelling.Schelling(args.size, args.size, args.empty_ratio,
      args.distribution, args.alpha, args.mean, args.std, args.max_iters,
      args.update_rate, args.num_races)
  model.populate()
  model.update()
  print 'Segregation index is ' + str(model.calculate_similarity())
  #model.plot('Grid')
  #plt.show()


if __name__ == "__main__":
  main()

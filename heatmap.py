#import matplotlib.pyplot as plt
import argparse
import schelling
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description="")
parser.add_argument("output_file", help="The location of output heatmap pdf")
parser.add_argument("-s" , "--size", dest="size", default=30, type=int,
                    help="size of the grid.")
parser.add_argument("-e" , "--empty_ratio", dest="empty_ratio", default=0.3, type=float,
                    help="Empty ratio of the grid")
parser.add_argument("-st", "--std", dest="std", default=0.1, type=float,
                    help="The std of guassian distribution")
parser.add_argument("-i", "--max_iterations", dest="max_iters", default=600, type=int,
                    help="max number of iterations")
parser.add_argument("-r", "--num_races", dest="num_races", default=2,
                    help="number of races")
parser.add_argument("-t", "--num_trials", dest="num_trials", default=2, type=int,
                    help="number of trials")
args = parser.parse_args()

def main():
  alpha = 0
  distribution = "Gaussian"
  means = np.arange(0.1, 1, 0.05)
  update_rates = np.arange(0, 1, 0.05)
  results = np.zeros((len(means) * len(update_rates),3))
  idx = 0
  for mean in means:
    for update_rate in update_rates:
      print 'Trying mean {} and update rate {}'.format(mean, update_rate)
      if update_rate == 0 and mean < 0.3:
        index = 1
      else:
        sum_index = 0.0
        for trial in range(0,args.num_trials):
          model = schelling.Schelling(args.size, args.size, args.empty_ratio,
                                      distribution, alpha, mean, args.std, args.max_iters,
                                      update_rate, args.num_races)
          model.populate()
          model.update()
          sum_index += model.calculate_similarity()
          index = sum_index/args.num_trials

      results[idx] = [mean, update_rate, index]
      idx += 1

  data = pd.DataFrame(results,
                      columns = ['Mean Tolerance' , 'Learning Rate', 'Segregation Index'])
  data = data.pivot(index='Mean Tolerance',
                    columns='Learning Rate', values='Segregation Index')
  fig, ax = plt.subplots(figsize=(12,12))
  sns.heatmap(ax = ax, data = data)
  sns.plt.show()
  fig.savefig(args.output_file)

if __name__ == "__main__":
  main()


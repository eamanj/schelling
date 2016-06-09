#import matplotlib.pyplot as plt
import itertools
import sys
import random
import copy
import argparse
import scipy.stats as stat

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


class Schelling:
  def __init__(self, width, height, empty_ratio, tolerance_distribution, alpha, mean, std,
               n_iterations, update_rate, races):
    self.width = width 
    self.height = height 
    self.races = races
    self.empty_ratio = empty_ratio
    self.tolerance_distribution = tolerance_distribution
    self.n_iterations = n_iterations
    self.empty_houses = []
    # agents is a dict from tuple (x_coord, y_coord) to list [race, tolerance]
    self.agents = {}
    self.alpha = alpha
    self.mean = mean
    self.std = std
    self.update_rate = update_rate

  def sample(self, tolerance_distribution, alpha, mean, std):
    if tolerance_distribution == 'Gaussian':
      return stat.norm.rvs(loc=mean,scale=std)
    elif tolerance_distribution == 'Powerlaw':
      return stat.powerlaw.rvs(alpha)
    else:
      sys.exit('Bad distribution ' + tolerance_distribution)

  def populate(self):
    self.all_houses = list(itertools.product(range(self.width),range(self.height)))
    random.shuffle(self.all_houses)

    self.n_empty = int( self.empty_ratio * len(self.all_houses) )
    self.empty_houses = self.all_houses[:self.n_empty]

    self.remaining_houses = self.all_houses[self.n_empty:]
    houses_by_race = [self.remaining_houses[i::self.races] for i in range(self.races)]

    for race_house in range(len(houses_by_race)):
      for address in houses_by_race[race_house]:
        tolerance = self.sample(self.tolerance_distribution, self.alpha, self.mean, self.std)
        self.agents[address] = [race_house + 1, tolerance]

  def count_similar_different_neighbors(self, agent_loc):
    x = agent_loc[0]
    y = agent_loc[1]
    race = self.agents[(x,y)][0]
    count_similar = 0
    count_different = 0

    if x > 0 and y > 0 and (x-1, y-1) not in self.empty_houses:
      if self.agents[(x-1, y-1)][0] == race:
        count_similar += 1
      else:
        count_different += 1
    if y > 0 and (x,y-1) not in self.empty_houses:
      if self.agents[(x,y-1)][0] == race:
        count_similar += 1
      else:
        count_different += 1
    if x < (self.width-1) and y > 0 and (x+1,y-1) not in self.empty_houses:
      if self.agents[(x+1,y-1)][0] == race:
        count_similar += 1
      else:
        count_different += 1
    if x > 0 and (x-1,y) not in self.empty_houses:
      if self.agents[(x-1,y)][0] == race:
        count_similar += 1
      else:
        count_different += 1        
    if x < (self.width-1) and (x+1,y) not in self.empty_houses:
      if self.agents[(x+1,y)][0] == race:
        count_similar += 1
      else:
        count_different += 1
    if x > 0 and y < (self.height-1) and (x-1,y+1) not in self.empty_houses:
      if self.agents[(x-1,y+1)][0] == race:
        count_similar += 1
      else:
        count_different += 1        
    if x > 0 and y < (self.height-1) and (x,y+1) not in self.empty_houses:
      if self.agents[(x,y+1)][0] == race:
        count_similar += 1
      else:
        count_different += 1        
    if x < (self.width-1) and y < (self.height-1) and (x+1,y+1) not in self.empty_houses:
      if self.agents[(x+1,y+1)][0] == race:
        count_similar += 1
      else:
        count_different += 1

    return (count_similar, count_different)

  def is_unsatisfied(self, agent_loc):
    tolerance = self.agents[agent_loc][1]
    # agent_loc is a list of form [x,y] corresponding to coordinates of the agent on the grid
    (count_similar, count_different) = self.count_similar_different_neighbors(agent_loc)
    if (count_similar+count_different) == 0:
      return False
    else:
      return float(count_different)/(count_similar + count_different) > tolerance
  
  def update_agent_tolerance(self, agent_loc):
    old_tolerance = self.agents[agent_loc][1]
    # should we update preferences?
    if not self.update_rate:
      return old_tolerance
    if self.update_rate > 1 or self.update_rate < 0:
      sys.exit('Bad update rate ' + self.update_rate)
    (count_similar, count_different) = self.count_similar_different_neighbors(agent_loc)
    if (count_similar+count_different) == 0:
      return old_tolerance
    different_neighbors_fraction = float(count_different)/(count_similar + count_different)
    new_tolerance = (1 - self.update_rate) * old_tolerance + self.update_rate * different_neighbors_fraction
    return new_tolerance

  def update(self):        
    for i in range(self.n_iterations):
      self.old_agents = copy.deepcopy(self.agents)
      n_changes = 0
      for agent_loc in self.old_agents:
        # compute new tolerance BEFORE MOVING, BUT UPDATE after moving. persumably, the
        # effect of interaction with neighbors shows up in the future, not in the current
        # time step
        new_agent_tolerance = self.update_agent_tolerance(agent_loc)
        if self.is_unsatisfied(agent_loc):
          agent_race = self.agents[agent_loc][0]
          empty_house = random.choice(self.empty_houses)
          # use new tolerance in new location
          self.agents[empty_house] = [agent_race, new_agent_tolerance]
          del self.agents[agent_loc]
          self.empty_houses.remove(empty_house)
          self.empty_houses.append(agent_loc)
          n_changes += 1
      if n_changes == 0:
        return
    # Out of the loop. it did not converge!
    print 'ERROR: grid did not converge in {} iterations'.format(self.n_iterations)
    print 'ERROR: increase either max number of iterations or mean tolerance of population'


  def move_to_empty(self, x, y):
    race = self.agents[(x,y)]
    empty_house = random.choice(self.empty_houses)
    self.updated_agents[empty_house] = race
    del self.updated_agents[(x, y)]
    self.empty_houses.remove(empty_house)
    self.empty_houses.append((x, y))

  def plot(self, title):
    fig, ax = plt.subplots()
    #If you want to run the simulation with more than 7 colors, you should set agent_colors accordingly
    agent_colors = {1:'b', 2:'r', 3:'g', 4:'c', 5:'m', 6:'y', 7:'k'}
    for agent in self.agents:
      ax.scatter(agent[0]+0.5, agent[1]+0.5, color=agent_colors[self.agents[agent][0]])
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_xlim([0, self.width])
    ax.set_ylim([0, self.height])
    ax.set_xticks([])
    ax.set_yticks([])
    #plt.savefig(file_name)

  def calculate_similarity(self):
    similarity = []
    for agent in self.agents:
      count_similar = 0
      count_different = 0
      x = agent[0]
      y = agent[1]
      race = self.agents[(x,y)][0]
      if x > 0 and y > 0 and (x-1, y-1) not in self.empty_houses:
        if self.agents[(x-1, y-1)][0] == race:
          count_similar += 1
        else:
          count_different += 1
      if y > 0 and (x,y-1) not in self.empty_houses:
        if self.agents[(x,y-1)][0] == race:
          count_similar += 1
        else:
          count_different += 1
      if x < (self.width-1) and y > 0 and (x+1,y-1) not in self.empty_houses:
        if self.agents[(x+1,y-1)][0] == race:
          count_similar += 1
        else:
          count_different += 1
      if x > 0 and (x-1,y) not in self.empty_houses:
        if self.agents[(x-1,y)][0] == race:
          count_similar += 1
        else:
          count_different += 1        
      if x < (self.width-1) and (x+1,y) not in self.empty_houses:
        if self.agents[(x+1,y)][0] == race:
          count_similar += 1
        else:
          count_different += 1
      if x > 0 and y < (self.height-1) and (x-1,y+1) not in self.empty_houses:
        if self.agents[(x-1,y+1)][0] == race:
          count_similar += 1
        else:
          count_different += 1        
      if x > 0 and y < (self.height-1) and (x,y+1) not in self.empty_houses:
        if self.agents[(x,y+1)][0] == race:
          count_similar += 1
        else:
          count_different += 1        
      if x < (self.width-1) and y < (self.height-1) and (x+1,y+1) not in self.empty_houses:
        if self.agents[(x+1,y+1)][0] == race:
          count_similar += 1
        else:
          count_different += 1
      try:
        similarity.append(float(count_similar)/(count_similar+count_different))
      except:
        similarity.append(1)
    return sum(similarity)/len(similarity)


def main():
  schelling = Schelling(args.size, args.size, args.empty_ratio,
      args.distribution, args.alpha, args.mean, args.std, args.max_iters,
      args.update_rate, args.num_races)
  schelling.populate()
  schelling.update()
  print 'Segregation index is ' + str(schelling.calculate_similarity())
  #schelling.plot('Grid')
  #plt.show()


if __name__ == "__main__":
  main()

import matplotlib.pyplot as plt
import itertools
import random
import copy
import scipy.stats as stat

#Hi
class Schelling:
	def __init__(self, width, height, empty_ratio, similarity_threshold, n_iterations, races = 2):
		self.width = width 
		self.height = height 
		self.races = races
		self.empty_ratio = empty_ratio
        #self.similarity_threshold = similarity_threshold
		#self.tolerance_distribution = 
		self.n_iterations = n_iterations
		self.empty_houses = []
		self.agents = {}

	def populate(self):
		self.all_houses = list(itertools.product(range(self.width),range(self.height)))
		random.shuffle(self.all_houses)

		self.n_empty = int( self.empty_ratio * len(self.all_houses) )
		self.empty_houses = self.all_houses[:self.n_empty]

		self.remaining_houses = self.all_houses[self.n_empty:]
		houses_by_race = [self.remaining_houses[i::self.races] for i in range(self.races)]
		
		for race_house in range(len(houses_by_race)):
			for address in houses_by_race[race_house]:
				self.agents[address] = [race_house + 1,stat.norm.rvs(loc=.8,scale=0.1)] # Second parameter is tolerance sampled from Gaussian
						
	def is_unsatisfied(self, x, y):
		race = self.agents[(x,y)][0]
		tolerance = self.agents[(x,y)][1]
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

		if (count_similar+count_different) == 0:
			return False
		else:
			return float(count_different)/(count_similar + count_different) > tolerance #self.similarity_threshold

	def update(self):        
		for i in range(self.n_iterations):
			self.old_agents = copy.deepcopy(self.agents)
			n_changes = 0
			for agent in self.old_agents:
				if self.is_unsatisfied(agent[0], agent[1]):
					agent_race = self.agents[agent]
					empty_house = random.choice(self.empty_houses)
					self.agents[empty_house] = agent_race
					del self.agents[agent]
					self.empty_houses.remove(empty_house)
					self.empty_houses.append(agent)
					n_changes += 1
			#print n_changes
			if n_changes == 0:
				break

	def move_to_empty(self, x, y):
		race = self.agents[(x,y)]
		empty_house = random.choice(self.empty_houses)
		self.updated_agents[empty_house] = race
		del self.updated_agents[(x, y)]
		self.empty_houses.remove(empty_house)
		self.empty_houses.append((x, y))
		
	def plot(self, title, file_name):
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
		
schelling = Schelling(30, 30, 0.1, 0.1, 500, 2)
schelling.populate()
schelling.plot('First plot','a')
plt.show()
schelling.update()
print 'segregation index: ' + str(schelling_1.calculate_similarity())
schelling.plot('Last plot', 'a')
plt.show()

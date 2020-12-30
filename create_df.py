import pandas as pd
import csv
from random import randint
#df = pd.read_csv("chess_year.csv", sep=' ')

keyDict = {"u", "v", "t", "w"}
temporal_dict = dict([(key, []) for key in keyDict])

with open('chess_year', newline='') as games:
    game_reader = csv.reader(games, delimiter='\t')
    for game in game_reader:
        temporal_dict['u'].append(game[0])
        temporal_dict['v'].append(game[1])
        temporal_dict['t'].append(game[2])
        temporal_dict['w'].append(randint(0,500))

df = pd.DataFrame(data=temporal_dict, columns=['u', 'v', 't','w'])
df.to_csv('weighted_chess', sep='\t', index=False)


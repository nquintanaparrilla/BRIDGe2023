import os
import matplotlib.pyplot as plt
import pandas as pd


folder_path = '/path/to/sc_mb_samples'
file_list = os.listdir(folder_path)
output_folder = '/path/to/sc_mb_samples'
output_filename = 'total_results.csv' 

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

data_frames = []

for i, file in enumerate(file_list):
    file_path = os.path.join(folder_path, file)
    if os.path.isfile(file_path):
        df = pd.read_csv(file_path) 
        data_frames.append(df)

concatenated_df = pd.concat(data_frames, ignore_index=True)

output_file_path = os.path.join(output_folder, output_filename)
concatenated_df.to_csv(output_file_path, index=False) 


s = {'S': 0}
u = {'U': 0}
a = {'A': 0}
threeT = {'3T': 0}
fiveT = {'5T': 0}
farthree = {'F3': 0}
farfive = {'F5': 0}

df = pd.read_csv('/path/to/total_results.csv')
for index, row in df.iterrows():
    reads = row['reads']
    status = row['status']
    
    if 'S' in status:
        s['S'] += reads
    elif 'U' in status:
        u['U'] += reads
    elif 'A' in status:
        a['A'] += reads
    elif '3T' in status:
        threeT['3T'] += reads
    elif '5T' in status:
        fiveT['5T'] += reads
    elif 'F3' in status:
        farthree['F3'] += reads
    elif 'F5' in status:
        farfive['F5'] += reads

definite = sum(s.values()) + sum(u.values()) + sum(a.values()) + sum(threeT.values()) + sum(fiveT.values()) 
total = definite + sum(farthree.values()) + sum(farfive.values())
sizes = [definite, farthree['F3'], farfive['F5']]
labels = ['definite categories', 'Far 3T', 'Far 5T']
colors = ['orchid', 'blue', 'green']
explode = (0.05, 0.05, 0.05)

fig1, ax1 = plt.subplots()
ax1.pie(sizes, colors=colors, labels=labels,autopct="%1.1f%%", startangle=90)

centre_circle = plt.Circle((0, 0), 0.70, fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)

ax1.legend(labels, title="Total Reads: {}".format(total), loc="center left", bbox_to_anchor=(1, 0.5))
ax1.axis('equal')
#plt.tight_layout()
#output_file_path = '/path/to/mb_sc_pie_chart1.png'
#plt.savefig(output_file_path)
plt.show()

s = {'S': 0}
u = {'U': 0}
a = {'A': 0}
threeT = {'3T': 0}
fiveT = {'5T': 0}

df = pd.read_csv('/path/to/sc_mb_samples/total_results.csv')
for index, row in df.iterrows():
    reads = row['reads']
    status = row['status']
    
    if 'S' in status:
        s['S'] += reads
    elif 'U' in status:
        u['U'] += reads
    elif 'A' in status:
        a['A'] += reads
    elif '3T' in status:
        threeT['3T'] += reads
    elif '5T' in status:
        fiveT['5T'] += reads

total = sum(s.values()) + sum(u.values()) + sum(a.values()) + sum(threeT.values()) + sum(fiveT.values())
sizes = [s['S'], u['U'], a['A'], threeT['3T'], fiveT['5T']]
labels = ['S', 'U', 'A', "3'UTR", "5'UTR"]
colors = ['gold', 'black', '#DF2945', 'cornflowerblue', 'yellowgreen']
explode = (0.05, 0.05, 0.05, 0.05, 0.05)

fig1, ax1 = plt.subplots()
ax1.pie(sizes, colors=colors, labels=labels, startangle=90)

centre_circle = plt.Circle((0, 0), 0.70, fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)
labels1 = ['Spliced {}'.format(s), 'Unspliced {}'.format(u), 'Ambiguous{}'.format(a), "3'UTR {}".format(threeT), "5'UTR {}".format(fiveT)]
ax1.legend(labels = labels1, title="Total Reads: {}".format(total), loc="center left", bbox_to_anchor=(1, 0.5))
ax1.axis('equal')
#plt.tight_layout()
#output_file_path = '/path/to/mb_sc_pie_chart2.png'
#plt.savefig(output_file_path)
plt.show()

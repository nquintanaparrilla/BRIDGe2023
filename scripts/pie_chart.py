import os
import pandas as pd
import matplotlib.pyplot as plt

# Folder path containing CSV files
folder_path = '/path/results/sc_mb_results'

# Initialize dictionaries to accumulate counts
s = {'S': 0}
u = {'U': 0}
a = {'A': 0}
i = {'I': 0}
threeT = {'3T': 0}
fiveT = {'5T': 0}

# Loop through each file in the folder
for filename in os.listdir(folder_path):
    if filename.endswith('.csv'):
        file_path = os.path.join(folder_path, filename)
        
        # Read CSV file
        df = pd.read_csv(file_path)
        print(f"Processing file: {file_path}")
        # Process each row in the DataFrame
        for index, row in df.iterrows():
            reads = row['reads']
            status = row['status']
            
            if 'S' in status:
                s['S'] += reads
            elif 'U' in status:
                u['U'] += reads
            elif 'A' in status:
                a['A'] += reads
            elif 'I' in status:
                i['I'] += reads
            elif '3T' in status:
                threeT['3T'] += reads
            elif '5T' in status:
                fiveT['5T'] += reads

# Calculate total reads and create the pie chart
total = sum(s.values()) + sum(u.values()) + sum(a.values()) + sum(i.values()) + sum(threeT.values()) + sum(fiveT.values())
sizes = [s['S'], u['U'], a['A'], i['I'], threeT['3T'], fiveT['5T']]
labels = ['S', 'U', 'Ambiguous within exon', 'Intronic', "3'UTR", "5'UTR"]
colors = ['gold', 'black', '#DF2945', 'orchid', 'cornflowerblue', 'yellowgreen']
explode = (0.05, 0.05, 0.05, 0.05, 0.05, 0.05)

fig1, ax1 = plt.subplots()
ax1.pie(sizes, colors=colors, autopct="%1.1f%%", labels=labels, startangle=90)

centre_circle = plt.Circle((0, 0), 0.70, fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)

ax1.legend(labels, title="Total Reads: {}".format(total), loc="lower left")
ax1.axis('equal')
plt.tight_layout()
output_file_path = '/path/mb_sc_pie_chart1.png'
plt.savefig(output_file_path)
plt.show()

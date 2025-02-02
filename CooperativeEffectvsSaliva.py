import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional Consensus file")
parser.add_argument("-o", "--output", help="Output folder")
args = parser.parse_args()

df = pd.read_csv(args.file, sep="\t")
# Create scatterplot with regression
sns.set(style='whitegrid')
plot = sns.lmplot(x='Saliva/SG', y='ratio_2', data=df, ci=None)

# Get the regression values
slope, intercept = plot.ax.get_lines()[0].get_data()
regression_equation = f'y = {slope:.2f}x + {intercept:.2f}'
print("Regression Equation:", regression_equation)

# Show the plot
plt.title('Scatterplot with Regression Line')
plt.savefig(args.output)
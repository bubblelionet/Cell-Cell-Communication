
# # This script fits a linear regression for each LR pair and combine to a big table
import numpy as np
import csv
import pickle
import matplotlib
import math
import pandas as pd
import matplotlib
from sklearn.utils import resample
from sklearn import linear_model
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from collections import Counter
from scipy.stats import chi2

def readCsv(x):
  """Parse file."""
  #colNames = ["method", "benchmark", "start", "end", "time", "memory"]
  df = pd.read_csv(x, sep=",")

  return df

def preprocessDf(df):
  """Transform ligand and receptor columns."""
  df["ligand-receptor"] = df["ligand"] + '-' + df["receptor"]
  df["component"] = df["component"] #.astype(str).str.zfill(2)

  return df


# Load subtype label
subtype_label_file=''
subtype_abundance_df = readCsv(subtype_label_file)
# subtype_label=[]
# with open(subtype_label_file) as file:
#     csv_file = csv.reader(file, delimiter=",")
#     for line in csv_file:
#         subtype_label.append(line)

# barcode_subtype=dict()
# for i in range(1,len(subtype_label)):
#     barcode_subtype[subtype_label[i][0]]= subtype_label[i][1]

# Load NEST output 
df = readCsv("/Users/victoriagao/local_docs/NEST/output/From_Fatema/exp2_B1_top20percent.csv")
# df = readCsv("/Users/victoriagao/local_docs/NEST/output/From_Fatema/exp1_C1_top20percent.csv")
output_processed = preprocessDf(df)


# ### Build feature matrix

### Merge NEST output with subtype label, and filter out the spots that are not in the subtype label
matched_spots_df = pd.merge(output_processed, subtype_abundance_df, left_on='from_cell', right_on='SpotID') # Change from_cell to to_cell if interested in the receptors


# filter out the LR that only appeared once
matched_spots_df = matched_spots_df[matched_spots_df['ligand-receptor'].duplicated(keep=False)] 
# Take only top 90% LR by frequency
lr_counts = matched_spots_df['ligand-receptor'].value_counts()
threshold = lr_counts.quantile(0.50)  # gives the value at the 50th percentile
top_percent_lrs = lr_counts[lr_counts >= threshold].index
matched_spots_df = matched_spots_df[matched_spots_df['ligand-receptor'].isin(top_percent_lrs)]
# Delete some columns
matched_spots_df = matched_spots_df.drop(columns=['to_cell', 'ligand', 'receptor', 'attention_score', 'component', 'from_id','to_id','SpotID'])


matched_spots_df

len(matched_spots_df['ligand-receptor'].unique())


# ### Fit logistic regression and outputs a big coefficient table for all regressions

unique_lr_pairs = matched_spots_df['ligand-receptor'].unique() # Get unique ligand-receptor pairs

results = []

# Iterate through each unique ligand-receptor pair
for lr_pair in unique_lr_pairs:
    # Prepare the feature matrix X and the target vector y
    # get all the columns from matched_spots_df except for the 'from_cell', 'edge_rank' and 'ligand-receptor'
    X_log_reg = matched_spots_df.drop(columns=["from_cell", "edge_rank", "ligand-receptor"])
    y_binary = ["yes" if lr == lr_pair else "no" for lr in matched_spots_df["ligand-receptor"]]
    # print dimensions of X and y
    # print(X_log_reg.shape, len(y_binary))

    # print(Counter(y_binary))
    
    # Build and fit the logistic model
    model_log_reg = linear_model.LogisticRegression(solver='lbfgs')
    model_log_reg.fit(X_log_reg, y_binary)
    
    # Extract coefficients and score
    coef = model_log_reg.coef_[0]  # Coefficients for the features
    score = model_log_reg.score(X_log_reg, y_binary)  # Accuracy score for how the model is fitted
    
    # Append the results (including LR pair, coefficients, and score) to our results list
    results.append([lr_pair] + list(coef) + [score])

# Define the column names for our results DataFrame
columns = ['Ligand-Receptor'] + X_log_reg.columns.tolist() + ['Accuracy Score']

# Create a DataFrame from our results
results_df = pd.DataFrame(results, columns=columns)
results_df

# Save the results to a CSV file
results_df.to_csv("/Users/victoriagao/local_docs/NEST/stored_variables/Celltype_LR_invidual_LogisticRegressions/logistic_regression_results.csv", index=False)



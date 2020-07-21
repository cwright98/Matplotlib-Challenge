#!/usr/bin/env python
# coding: utf-8

# ## Observations and Insights 

# 

# In[2]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
df=mouse_metadata.merge(study_results, on="Mouse ID", how="outer")

# Display the data table for preview

print(len(df.index))


# In[3]:


# Checking the number of mice.
df['Mouse ID'].nunique()


# In[4]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint.
df=df.drop_duplicates(subset=[ "Mouse ID", "Timepoint"])


# In[5]:


# Optional: Get all the data for the duplicate mouse ID. 


# In[6]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.


# In[7]:


# Checking the number of mice in the clean DataFrame.
df['Mouse ID'].nunique()


# ## Summary Statistics

# In[8]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# This method is the most straighforward, creating multiple series and putting them all together at the end.


df.groupby('Drug Regimen').agg({'Tumor Volume (mm3)': ['mean', 'median', 'var', 'std']})


# In[9]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# This method produces everything in a single groupby function
#are u asking me to do the same thing twice? I'm confused if this is the same thing as above or not? It sounds like the same thing


# ## Bar and Pie Charts

# In[66]:


# Generate a bar plot showing the total number of mice for each treatment throughout the course of the study using pandas. 


# In[63]:


# Generate a bar plot showing the total number of mice for each treatment throughout the course of the study using pyplot.

drug_count=df.groupby('Drug Regimen').count()
drug_count=drug_count["Mouse ID"]
x_axis=drug_count.index

plt.bar(x_axis, drug_count)
plt.xticks(rotation=45)


# In[12]:


# Generate a pie plot showing the distribution of female versus male mice using pandas



# In[13]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
labels=["Female", "Male"]
sizes=df.groupby('Sex').count()
sizes['Mouse ID']
graph_sizes=[sizes['Mouse ID'][0], sizes['Mouse ID'][1]]
graph_sizes
plt.pie(graph_sizes, labels=labels)


# ## Quartiles, Outliers and Boxplots

# In[17]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin


# Start by getting the last (greatest) timepoint for each mouse
#find the max timepoint here, it's not always 45
last_timepoint=df.groupby('Mouse ID')
last_timepoint=last_timepoint['Timepoint'].max()
new_df=df.groupby(['Mouse ID', 'Drug Regimen'])['Timepoint'].max().to_frame()
# Merge this group df with the original dataframe to get the tumor volume at the last timepoint
step2=new_df.merge(df, on=['Timepoint', 'Mouse ID'])
step3=step2.loc[(step2['Drug Regimen']=='Capomulin') | (step2['Drug Regimen']=='Ramicane') | (step2['Drug Regimen']=='Infubinol')
               | (step2['Drug Regimen']=='Ceftamin')]
step4=step3[['Timepoint', 'Mouse ID', 'Drug Regimen', 'Tumor Volume (mm3)']]
step4=step4.sort_values('Drug Regimen')
step4.head()


# In[18]:


# Put treatments into a list for for loop (and later for plot labels)
drug_names=["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Create empty list to fill with tumor vol data (for plotting)
Capomulin=step4.loc[step4['Drug Regimen']=='Capomulin']['Tumor Volume (mm3)']
Ramicane=step4.loc[step4['Drug Regimen']=='Ramicane']['Tumor Volume (mm3)']
Infubinol=step4.loc[step4['Drug Regimen']=='Infubinol']['Tumor Volume (mm3)']
Ceftamin=step4.loc[step4['Drug Regimen']=='Ceftamin']['Tumor Volume (mm3)']
drug_list=[Capomulin, Ramicane, Infubinol, Ceftamin]
fig1, ax1 = plt.subplots()
#plt.xticks(np.arange(1, 5), drug_names)

ax1.set_title('Final Tumor Size for top 4 drugs')
ax1.set_ylabel('Size (mm3)')
ax1.set_xlabel("Drug regmimen")
ax1.boxplot([Capomulin, Ramicane, Infubinol, Ceftamin])

plt.show()
for drug_index in range(4):
    quartiles = drug_list[drug_index].quantile([.25,.5,.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq-lowerq
    print(f"Here are the results for {drug_names[drug_index]}!")

    print(f"The lower quartile of tumor size is: {lowerq}")
    print(f"The upper quartile of tumor size is: {upperq}")
    print(f"The interquartile range of tumor size is: {iqr}")
    print(f"The the median of tumor size is: {quartiles[0.5]} ")

    lower_bound = lowerq - (1.5*iqr)
    upper_bound = upperq + (1.5*iqr)
    print(f"Values below {lower_bound} could be outliers.")
    print(f"Values above {upper_bound} could be outliers.")


# In[19]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest


# ## Line and Scatter Plots

# In[36]:


# Generate a line plot of time point versus tumor volume for a mouse treated with Capomulin
Cap=df[['Tumor Volume (mm3)', 'Timepoint', 'Drug Regimen']].loc[df['Drug Regimen']=='Capomulin']
Cap_var=Cap.groupby('Timepoint').mean()['Tumor Volume (mm3)']

fig, ax=plt.subplots(figsize=(20, 6))
ax.set_ylabel('Tumor Size')
ax.set_title('Timepoint vs tumor volume')
Cap_var.plot()

plt.show()


# In[ ]:





# In[37]:


# Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin regimen
CapScatter=df.loc[df['Drug Regimen']=='Capomulin']
CapScatter.groupby('Weight (g)')['Tumor Volume (mm3)'].mean().plot()


# ## Correlation and Regression

# In[58]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen
df.head()
CapScatter=CapScatter[['Tumor Volume (mm3)', 'Weight (g)']]
x_values=CapScatter['Weight (g)']
y_values=CapScatter['Tumor Volume (mm3)']


# In[59]:


plt.scatter(x_values, y_values)


# In[62]:


from scipy.stats import linregress
(slope, intercept, rvalue, pvalue, stderr) = linregress(x_values, y_values)
regress_values = x_values * slope + intercept
line_eq = "y = " + str(round(slope,2)) + "x + " + str(round(intercept,2))
plt.scatter(x_values,y_values)
plt.plot(x_values,regress_values,"r-")
plt.annotate(line_eq,(6,10),fontsize=15,color="red")
plt.xlabel('')
plt.ylabel('Tumor Size')
plt.show()


# In[ ]:





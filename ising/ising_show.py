import pandas as pd
import matplotlib.pyplot as plt


# read the CSV file using pandas.read_csv()
df = pd.read_csv('data.csv', header=[0, 1])

# extract the arrays from the DataFrame
temps = df["array1"].values
abs_mag = df["array2"].values

plt.plot(temps, abs_mag,"*")
plt.show()

import matplotlib.pyplot as plt

# Define the file path
file_path = "/home/manoj/vins/squeakr/data/filter_data_10k.txt"

# Initialize an empty list to store the counts
counts = []

# Read the file line by line
with open(file_path, "r") as f:
    for line in f:
        # Split the line into the string and count
        string, count = line.strip().split("\t")
        # Convert the count to an integer and append to the counts list
        counts.append(int(count))

# Plot a histogram of the counts
plt.hist(counts, bins=50)
plt.xlabel("Count")
plt.ylabel("Number of Strings")
# Save the plot as a PNG file
plt.savefig("hist-output.png")
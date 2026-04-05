"""
Name: A.D.D.C. Wijethunge
Index: s16562
PPI Network Neighborhood Analyzer
A program to analyze protein- protein interaction networks
Time duration: 15/01/ 2026 to 01/04/ 2026
"""

import networkx as nx #dreate and analyze networks
import pandas as pd # for reading data file
import numpy as np # for math operation
import os # for file operation
import matplotlib.pyplot as plt
import seaborn as sns
import re

# Set seaborn style for all plots
sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1.2)

class PPIAnalyzer:

    #create objects
    def __init__(self):
        self.graph=None # store network
        self.annotations={} # store protein functions
        self.function_proteins=[] # list of known functional proteins
        self.function_categories={} #function categories if available
        self.best_n= 1 #best neighborhood distance(starting default value)



# ============================================================
#   Load a PPI network file and creates a NetworkX graph
#=============================================================

    def load_ppi_network(self, file_path):
        try:
            # read the ppi file
            print("\n" + "=" * 60)
            print("📊 LOADING NETWORK")
            print("=" * 60)
            print("  📁 File: " + os.path.basename(file_path))

            #check whether the file exists
            if not os.path.exists(file_path):
                print("  ❌ File not found at " + file_path)
                return None

            #read files (.csv, .tsv, space separated files are supported)
            try:
                data = pd.read_csv(file_path, sep=r"\s+", engine="python")#space seperator is the file format(add 'r' for raw string)
            except:
                try:
                    data = pd.read_csv(file_path, sep="\t")
                except:
                    data = pd.read_csv(file_path) #default(comma)

            print("  ✅ Loaded {:,} interactions".format(len(data)))
            print("  📄 Data format: {} columns".format(len(data.columns)))

            #check at least have 2 columns
            if len(data.columns)<2:
                print("  ❌ File needs at least 2 columns(protein1, protein2)")
                return None

            #create an empty graph
            self.graph=nx.Graph()
            #edges_added=0

            #add each iteration as an edge
            for index, row in data.iterrows(): # use pandas
                try:
                    protein1= str(row.iloc[0])
                    protein2 = str(row.iloc[1])

                    #skip if either protein is NaN (missing values)
                    if pd.isna(protein1) or pd.isna(protein2):
                        continue

                    #add the connection
                    self.graph.add_edge(protein1.strip(), protein2.strip())
                    #edges_added +=1

                except Exception as e:
                    print("  ⚠️ Could not process row {}: {}".format(index, e))
                    continue

            print("\n  ✅ NETWORK LOADED SUCCESSFULLY")
            print("  🧬 Proteins: {:,}".format(self.graph.number_of_nodes()))
            print("  🔗 Interactions: {:,}".format(self.graph.number_of_edges()))

            return self.graph

        except Exception as e:
            print("  ❌ Could not load file: " + str(e))
            return None


# ============================================================
#   Network Statistics
#   (Count proteins and interactions)
# =============================================================

    def get_network_stats(self):
        #check if graph exists
        if self. graph is None:
            print("  ❌ No graph loaded! Load the network first.")
            return 0, 0 # return 2 int values as 0, when the graph isn't present

        #get basic counts
        num_proteins= self.graph.number_of_nodes() # how many proteins?
        num_interactions = self.graph.number_of_edges() # how many connections?

        #calculate average degree
        if num_proteins > 0:
            avg_degree = (2 * num_interactions) / num_proteins # average connection per protein
        else:
            avg_degree = 0

        #display results
        print("\n" + "=" * 60)
        print("📈 NETWORK SUMMARY")
        print("=" * 60)
        print("  🧬 Proteins              : {:,}".format(num_proteins))
        print("  🔗 Interactions          : {:,}".format(num_interactions))
        print("  📊 Average connections   : {:.1f}".format(avg_degree))
        print("=" * 60)
        return num_proteins, num_interactions



# ============================================================
#   Calculate Protein degree
#   (Count how many neighbors a protein has)
# =============================================================
    def calculate_degree(self, protein):
        if self.graph is None:
            print("  ❌ Graph is not Loaded!")
            return 0

        if protein in self. graph:
            degree= self.graph.degree(protein) # how many proteins does this one interact with?
            print("  🔗 Degree of protein {}: {}".format(protein, degree))
            return degree
        else:
            print("  ❌ Protein {} not found in the network!".format(protein))
            return 0


    #=================Degree distribution as a Histogram===================
    def plot_degree_distribution(self, save_path= None):
        if self.graph is None:
            print("  ❌ Graph is not Loaded!")
            return

        #get all degrees of each protein node
        degrees = [degree for node, degree in self.graph.degree()]

        print("  📊 Collected degree for {} proteins".format(len(degrees)))

        #create the plot
        plt.figure(figsize=(12,6))
        #create a histogram- matplotlib package
        plt.hist(degrees, bins=15, color='lightblue',edgecolor='black', alpha=0.7)
        #add labels and title
        plt.xlabel('Degree(Number of connections)',fontsize=12)
        plt.ylabel('Frequency(Number of Proteins)', fontsize=12)
        plt.title('Degree distribution of PPI Network', fontsize=14)
        plt.xticks(range(0, max(degrees) + 10, 10))
        plt.grid(True, alpha=0.3)

        #add statistics text
        stats_text = 'Total Proteins: {}\n'.format(len(degrees))
        stats_text += 'Mean Degree: {:.2f}\n'.format(np.mean(degrees))
        stats_text += 'Max Degree: {}\n'.format(max(degrees))
        stats_text += 'Min Degree: {}'.format(min(degrees))

        plt.text(0.95, 0.95, stats_text,
                 transform=plt.gca().transAxes,
                 verticalalignment= 'top',
                 horizontalalignment='right',
                 fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                 )

        #save/ show
        if save_path:
            #create directory if it doesn't exist
            directory = os.path.dirname(save_path)
            if directory:  # Only create if directory is not empty
                os.makedirs(directory, exist_ok= True)
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print("  💾 Plot saved to: {}".format(save_path))



# ============================================================
#   Load function proteins
#   (Read the file containing proteins that already have the function)
# =============================================================

    def load_function_proteins(self, file_path):
        """Load function proteins(seeds) from tsv file"""
        try:
            self.function_proteins = []

            print("\n" + "=" * 60)
            print("📂 LOADING FUNCTION PROTEINS")
            print("=" * 60)

            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split('\t')

                    # Extract protein IDs from all columns
                    all_ids = []
                    for part in parts:
                        # Split by common separators
                        for sep in ['|', ';', ',', ' ']:
                            if sep in part:
                                all_ids.extend(part.split(sep))
                            else:
                                all_ids.append(part)

                    # Clean and filter IDs
                    for pid in all_ids:
                        pid = pid.strip()

                        # Remove taxon prefixes (e.g., 4932.)
                        if '.' in pid and pid.split('.')[0].isdigit():
                            pid = pid.split('.')[-1]

                        # Remove version suffixes
                        if '.' in pid and pid.split('.')[-1].isdigit():
                            pid = pid.split('.')[0]

                        # Keep if it looks like a protein ID
                        if pid and len(pid) >= 3 and re.search(r'[A-Za-z]', pid):
                            self.function_proteins.append(pid)

            # Remove duplicates
            self.function_proteins = list(set(self.function_proteins))

            print("  ✅ Loaded {:,} unique function proteins".format(len(self.function_proteins)))
            return self.function_proteins

        except Exception as e:
            print("  ❌ Error loading function proteins: {}".format(e))
            return []


# ========================= NEIGHBORHOOD ANALYSIS METHODS ================================

# ============================================================
#   Count Annotated neighbors
#   (Count how many neighbors of a protein are functional)
# =============================================================

    def count_annotated_neighbors(self, protein):
        if self.graph is None:
            print("  ❌ Graph is not Loaded!")
            return 0

        if not self.function_proteins:
            print("  ❌ No function proteins loaded!")
            return 0

        if protein not in self.graph:
            print("  ❌ Protein {} not in network! ".format(protein))
            return 0

        #get all neighbors
        neighbors= list(self.graph.neighbors(protein))

        #convert function list to set
        func_set=set(self.function_proteins)

        #count neighbors in function set
        count=0
        annotated_neighbors=[]
        for neighbor in neighbors:
            if neighbor in func_set: #if this neighbor functional?
                count +=1
                annotated_neighbors.append(neighbor)

        #output results
        print("\n" + "=" * 60)
        print("🔍 NEIGHBOR ANALYSIS: {}".format(protein))
        print("=" * 60)
        print("  🔗 Total neighbors              : {}".format(len(neighbors)))
        print("  ✅ Annotated neighbors count    : {} ({:.1%})".format(count, count/len(neighbors)))
        if annotated_neighbors:
            # Show first 10 to avoid clutter
            display_neighbors = annotated_neighbors[:15]
            neighbor_str = ', '.join(display_neighbors)
            print("  📌 Annotated neighbors found   : {}".format(neighbor_str))
            if len(annotated_neighbors) > 15:
                print("  📌 ... and {} more".format(len(annotated_neighbors) - 15))

        # Recalculate p (same as Hishigaki)
        total_proteins = self.graph.number_of_nodes()
        func_set = set(self.function_proteins)
        functional = [p for p in self.graph.nodes() if p in func_set]
        p = len(functional) / total_proteins if total_proteins > 0 else 0
        expected = len(neighbors) * p
        print("  📊 Expected by chance           : {:.1f} ({:.1%})".format(expected, p))
        if count < expected:
            print("  📉 This is {:.1f} FEWER than expected".format(expected - count)) # it might have a specialized role or regulatory protein
        else:
            print("  📈 This is {:.1f} MORE than expected".format(count - expected))
        print("=" * 60)
        return count



# ============================================================
#   Getting n-neighbors
#   (Get all proteins within n steps of a protein)
# =============================================================

    def get_n_neighbors(self, protein, n=1):
        # get all proteins reachable within n interactions
        #basic check
        if self.graph is None or protein not in self.graph:
            return []

        # direct neighbors(n=1)
        if n == 1:
            return list(self.graph.neighbors(protein))

        #indirect connection (n>1)
        # use Breadth first search algorithm
        visited = set([protein])
        current_level = set([protein])

        for distance in range(n):
            next_level = set()
            for node in current_level:
                for neighbor in self.graph.neighbors(node):
                    if neighbor not in visited:
                        visited.add(neighbor)
                        next_level.add(neighbor)
            current_level = next_level

        #remove query protein itself
        visited.remove(protein)
        return list(visited)


    #=================Plot Annotation Distribution==================
    def plot_annotation_distribution(self, remove_zeros= False, save_path= None):
        if self.graph is None:
            print("  ❌ Graph is not Loaded!")
            return

        if not self.function_proteins:
            print("  ❌ No function proteins loaded!")
            return


        #calculate counts for all proteins
        counts=[]
        func_set= set(self.function_proteins)

        for protein in self.graph.nodes():
            neighbors= list(self.graph.neighbors(protein))

            #count how many neighbors in function list
            count=0
            for neighbor in neighbors:
                if neighbor in func_set:
                    count += 1
            #apply zero removal if required
            if not(remove_zeros and count==0):
                counts.append(count)

        print("  📊 Calculated counts for {} proteins".format(len(counts)))

        #plotting
        plt.figure(figsize=(12,6))
        #create a histogram
        plt.subplot(1,2,1)
        plt.hist(counts, bins=15, color="Lightgreen", edgecolor='black', alpha=0.7)
        plt.xlabel("Number of Annotated Neighbors", fontsize=12)
        plt.ylabel("Frequency", fontsize=12)
        plt.title("Histogram of Distribution of Annotated proteins in Neighborhoods", fontsize=14)
        plt.grid(True, alpha=0.3)

        #Set x-axis based on actual data instead of removing ticks
        max_count = max(counts)
        #Show ticks at reasonable intervals based on data range
        if max_count <= 10:
            step = 1
        elif max_count <= 20:
            step = 2
        elif max_count <= 50:
            step = 5
        else:
            step = 10

        ticks = range(0, max_count + step, step)
        plt.xticks(ticks)  #actual ticks
        plt.xlim(-0.5, max_count + 0.5)  #padding

        #main title
        title= "Distribution of Annotated Proteins in Neighborhoods"
        if remove_zeros:
            title +="(Zeros Removed)"
            plt.suptitle(title, fontsize=14, fontweight="bold")

        #add statistics
        stats_text="Total proteins: {}\n".format(len(counts))
        stats_text += "Mean count: {:.2f}\n".format(np.mean(counts))
        stats_text += "Max count: {}".format(max(counts))

        plt.figtext(0.95, 0.95, stats_text,
                    verticalalignment='top',
                    horizontalalignment='right',
                    fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                    )

        # save file
        if save_path:
            #create directory if it doesn't exist
            directory = os.path.dirname(save_path)
            if directory:  # Only create if directory is not empty
                os.makedirs(directory, exist_ok= True)
            plt.savefig(save_path, dpi=300, bbox_inches= 'tight')
            print("  💾 Plot saved to: {}".format(save_path))

        plt.close()  # Close the figure to free memory


# ============================================================
#   Plot degree distribution
#   (histogram shows how many proteins have each degree)
# =============================================================
    def create_degree_plot_figure(self):
        """Create and return a degree distribution figure for GUI"""
        if self.graph is None:
            print("  ❌ Graph is not Loaded!")
            return None

        # Create figure without showing
        fig, ax = plt.subplots(figsize=(10, 6))

        # Get all degrees
        degrees = [degree for node, degree in self.graph.degree()]

        # Create histogram
        ax.hist(degrees, bins=30, color='lightblue', edgecolor='black', alpha=0.7)
        ax.set_xlabel('Degree(Number of connections)', fontsize=12)
        ax.set_ylabel('Frequency(Number of Proteins)', fontsize=12) #how many proteins have that many connection
        ax.set_title('Degree distribution of PPI Network', fontsize=14)
        ax.set_xticks(range(0, max(degrees) + 10, 10))
        ax.grid(True, alpha=0.3)

        # Add statistics box (wheat color)
        stats_text = 'Total Proteins: {}\n'.format(len(degrees))
        stats_text += 'Mean Degree: {:.2f}\n'.format(np.mean(degrees))
        stats_text += 'Max Degree: {}\n'.format(max(degrees))
        stats_text += 'Min Degree: {}'.format(min(degrees))
        ax.text(0.95, 0.95, stats_text,
                transform=ax.transAxes,
                verticalalignment='top',
                horizontalalignment='right',
                fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        return fig



# ============================================================
#   Plot annotation distribution
#   (show how many annotated neighbors each protein has)
# =============================================================

    def create_annotation_plot_figure(self, remove_zeros=True):
        """Create and return annotation distribution figure for GUI - with boxplot"""
        if self.graph is None:
            print("  ❌ Graph is not Loaded!")
            return None

        if not self.function_proteins:
            print("  ❌ No function proteins loaded!")
            return None

        # Calculate counts
        func_set = set(self.function_proteins)
        counts = [] #store how many functional neighbors each protein has

        for protein in self.graph.nodes(): # look at every protein in network and get its friends
            neighbors = list(self.graph.neighbors(protein))
            count = 0

            for neighbor in neighbors:
                if neighbor in func_set: #count how many friends are functional
                    count += 1

            #decide whether to keep this count
            if not(remove_zeros and count == 0): #exclude proteins with 0 functional neighbors
                counts.append(count) # store count for each protein

        print("  📊 Calculated counts for {} proteins".format(len(counts)))


        # Create figure with 2 subplots (HISTOGRAM + BOXPLOT)
        fig = plt.figure(figsize=(12, 6))

        # Left panel: Histogram (shows distribution)
        ax1 = plt.subplot(1, 2, 1)
        ax1.hist(counts, bins=15, color="lightgreen", edgecolor='black', alpha=0.7)
        ax1.set_xlabel("Number of Annotated Neighbors", fontsize=12) # number of functional neighbors
        max_count = max(counts)
        ax1.set_xticks(range(0, max_count + 10, 10))  # ticks at 0, 10, 20, ..., 70
        ax1.set_xlim(-2, max_count + 5)  # paddings
        ax1.set_ylabel("Frequency", fontsize=12) # how many proteins have that many functional neighbors
        ax1.set_title("Histogram of Annotated Neighbors", fontsize=14)
        ax1.grid(True, alpha=0.3)

        # Right panel: Boxplot (show summary statistics)
        ax2 = plt.subplot(1, 2, 2)
        ax2.boxplot(counts, vert=True, patch_artist=True,
                    boxprops=dict(facecolor='lightgreen', alpha=0.7))
        ax2.set_ylabel('Annotated Neighbors Count', fontsize=12)
        ax2.set_title('Boxplot of Annotation Distribution', fontsize=14)
        ax2.grid(True, alpha=0.3)

        # Add statistics text
        stats_text = "Total proteins: {}\n".format(len(counts))
        stats_text += "Mean count: {:.2f}\n".format(np.mean(counts))
        stats_text += "Max count: {}".format(max(counts))
        plt.figtext(0.95, 0.95, stats_text,
                    verticalalignment='top',
                    horizontalalignment='right',
                    fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))


        # Add main title
        main_title = "Distribution of Annotated Proteins in Neighborhoods"
        if remove_zeros:
            main_title += " (Zeros Removed)"
        plt.suptitle(main_title, fontsize=14, fontweight='bold')
        return fig



# ============================================================
#   Hishigaki Algorithm score
#   (predict what protein might have the function)
# =============================================================

    # Predict candidate genes
    def calculate_hishigaki_scores(self):
        if self.graph is None:
            print("  ❌ Graph is not Loaded!")
            return {}

        if not self.function_proteins:
            print("  ❌ No function proteins loaded!")
            return {}

        #convert to set for faster lookup
        func_set= set(self.function_proteins) #seed protein set
        scores={} #empty dictionary to store scores {protein: score}

        #calculate background probability
        total_proteins = self.graph.number_of_nodes() #total proteins in network
        # functional proteins in the network
        total_functional = len([p for p in self.graph.nodes() if p in func_set])

        # probability a random protein has the function (between 0-1)
        if total_proteins > 0:
            p = total_functional / total_proteins
        else:
            p = 0

        print("\n" + "=" * 60)
        print("📊 CALCULATING χ² HISHIGAKI SCORES")
        print("=" * 60)
        print("  📈 Background probability          : p = {:.4f}".format(p))
        print("  🧬 Total proteins                  : {}".format(total_proteins))
        print("  📚 Annotated proteins              : {}".format(len(func_set)))
        print("  ✅ Functional proteins in network  : {}".format(total_functional))
        print("=" * 60)


        #calculate scores for each protein
        for protein in self.graph.nodes():
            if protein in func_set:  # skip already functional proteins
                scores[protein]= 0.0 # already known, not a candidate
                continue

            #get neighbors within best_n distance
            neighbors= self.get_n_neighbors(protein, self.best_n)

            #skip if no neighbor
            if not neighbors:
                scores[protein]= 0.0
                continue

            #count annotated neighbors
            annotated_count=0
            for neighbor in neighbors:
                if neighbor in func_set: # is this neighbor functional?
                    annotated_count+=1

            #calculated expected number of annotated neighbor
            total_neighbors = len(neighbors)
            expected = total_neighbors * p

            # Calculate χ² Hishigaki score (high score--> strong candidate)
            if expected > 0:
                # χ² = (observed - expected)² / expected
                chi_square = ((annotated_count - expected) ** 2) / expected
                scores[protein] = chi_square
            else:
                scores[protein] = 0.0

        print("\n  ✅ Calculated χ² Hishigaki scores for {:,} proteins".format(len(scores)))
        print("  🔗 Using neighborhood distance: n={}".format(self.best_n))
        print("  📊 Score range: [{:.2f}, {:.2f}]".format(min(scores.values()), max(scores.values())))
        print("-" * 60)

        return scores


# ============================================================
#   Plot Hishigaki distribution
#   (crate histogram of all scores)
# =============================================================

    def create_hishigaki_plot_figure(self, scores, remove_zeros=False):
        if not scores:
            print("  ❌ No scores calculated!")
            return None

        # Get score values
        score_values = list(scores.values())
        # Remove zeros if needed
        if remove_zeros:
            score_values = [s for s in score_values if s > 0]

        if not score_values:
            print("  ❌ No scores to plot!")
            return None

        # Create figure
        fig = plt.figure(figsize=(14, 6))

        # ----------LEFT: Histogram------------------
        ax1 = plt.subplot(1, 2, 1)
        ax1.hist(score_values, bins=10, color='salmon', edgecolor='black', alpha=0.7)
        ax1.set_xlabel('Hishigaki Score', fontsize=12)
        ax1.set_ylabel('Frequency', fontsize=12)
        ax1.set_title('Histogram of Hishigaki scores', fontsize=14)
        ax1.grid(True, alpha=0.3)
        # set x_axis limits and let matplotlib to choose ticks
        max_score = max(score_values)
        ax1.set_xlim(-0.5, max_score + 1)

        # Add mean and median lines
        mean_score = np.mean(score_values)
        median_score = np.median(score_values)
        ax1.axvline(mean_score, color='red', linestyle='--', linewidth=2, label='Mean: {:.2f}'.format(mean_score))
        ax1.axvline(median_score, color='green', linestyle='--', linewidth=2, label='Median: {:.2f}'.format(median_score))
        ax1.legend(fontsize=10, facecolor='wheat')



        #--------------RIGHT: Boxplot-------------------
        ax2 = plt.subplot(1, 2, 2)
        ax2.boxplot(score_values, vert=True, patch_artist=True,
                    boxprops=dict(facecolor='salmon', alpha=0.7))
        ax2.set_ylabel('Hishigaki Score', fontsize=12)
        ax2.set_title('Boxplot of Hishigaki scores', fontsize=14)
        ax2.grid(True, alpha=0.3)
        ax2.set_xticks([])  # Remove x-axis numbers


        plt.suptitle("Distribution of χ² Hishigaki scores (n={},  p={:.4f})".format(self.best_n, mean_score/100), fontsize=14, fontweight='bold')

        return fig



    #--------------------------------------------------------
    # Plotting Hishigaki score Distribution (Histogram Only)
    #--------------------------------------------------------
    def plot_hishigaki_distribution(self, scores, remove_zeros=False, save_path=None):
        """Display Hishigaki score distribution as a histogram only"""
        # Check for the hishgaki scores
        if not scores:
            print("  ❌ No scores calculated!")
            return

        # Get score values, remove zeros if needed
        score_values = list(scores.values())
        if remove_zeros:
            score_values = [s for s in score_values if s > 0]

        if not score_values:
            return

        # Create plot - histogram only
        plt.figure(figsize=(10, 6))
        plt.hist(score_values, bins=10, color='salmon', edgecolor='black', alpha=0.7)
        plt.xlabel('χ² Hishigaki Score', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
        plt.title('Histogram of χ² Hishigaki scores (n={})'.format(self.best_n), fontsize=14)
        plt.grid(True, alpha=0.3)
        # Set x-axis based on actual data
        max_score = max(score_values)
        # Let matplotlib choose ticks based on data range
        plt.xlim(-0.5, max_score + 1)


        # Add mean and median lines
        mean_score = np.mean(score_values)
        median_score = np.median(score_values)
        plt.axvline(mean_score, color='red', linestyle='--', label='Mean: {:.2f}'.format(mean_score))
        plt.axvline(median_score, color='green', linestyle='--', label='Median: {:.2f}'.format(median_score))
        plt.legend(fontsize=10)

        # Add statistics box
        stats_text = "Total: {}\nMean: {:.2f}\nMedian: {:.2f}\nMax: {:.2f}".format(len(score_values), mean_score, median_score, max(score_values))
        plt.text(0.95, 0.95, stats_text, transform=plt.gca().transAxes,
                 verticalalignment='top', horizontalalignment='right', fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        plt.tight_layout()

        #save if path given
        if save_path:
            folder = os.path.dirname(save_path)
            if folder:
                os.makedirs(folder, exist_ok=True)
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print("  💾 Histogram saved to: {}".format(save_path))
        plt.close()

    #-----------------------------------------------------
    # Plotting Hishigaki score Boxplot (Separate Method)
    #-----------------------------------------------------
    def plot_hishigaki_boxplot(self, scores, remove_zeros=False, save_path=None):
        """Display Hishigaki score distribution as a boxplot only"""
        if not scores:
            print("  ❌ No scores calculated!")
            return

        # Get scores, remove zeros if needed
        score_values = list(scores.values())
        if remove_zeros:
            score_values = [s for s in score_values if s > 0]

        if not score_values:
            return

        # Create plot - boxplot only
        plt.figure(figsize=(8, 6))
        plt.boxplot(score_values, vert=True, patch_artist=True,
                    boxprops=dict(facecolor='salmon', alpha=0.7))
        plt.ylabel('χ² Hishigaki Score', fontsize=12)
        # Add title
        title = 'Boxplot of χ² Hishigaki scores'
        if remove_zeros:
            title += ' (Non-zero only)'
        plt.title(title, fontsize=14)

        plt.grid(True, alpha=0.3)
        plt.xticks([])  # Remove x-axis numbers for single boxplot

        # Add statistics box
        stats_text = "Total: {}\nMean: {:.2f}\nMedian: {:.2f}\nMax: {:.2f}".format(len(score_values), np.mean(score_values), np.median(score_values), max(score_values))
        plt.text(0.95, 0.95, stats_text, transform=plt.gca().transAxes,
                 verticalalignment='top', horizontalalignment='right', fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        plt.tight_layout()

        # Save if path given
        if save_path:
            folder = os.path.dirname(save_path)
            if folder:
                os.makedirs(folder, exist_ok=True)
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print("  💾 Boxplot saved to: {}".format(save_path))
        plt.close()



# ============================================================
#   Save results to .csv file
#   (save scores to a .csv file)
# =============================================================

    def save_results(self, scores, save_path):
        if not scores:
            print("  ❌ No scores to save!")
            return

        try:
            #convert dictionary to data frame
            df=pd.DataFrame(list(scores.items()), columns=['Proteins','Hishigaki_Scores'])

            #sort by score(descending order)
            df=df.sort_values('Hishigaki_Scores', ascending =False) # best candidate is the top

            #create directory if needed
            directory = os.path.dirname(save_path)
            if directory:  # Only create if directory is not empty
                os.makedirs(directory, exist_ok= True)

            #save to csv(don't add row numbers to the file)
            df.to_csv(save_path, index=False)

            print("  💾 Results saved to: {}".format(save_path))
            print("  📊 Total records: {}".format(len(df)))

        except Exception as e:
            print("  ❌ Error: {}".format(e))



# ============================================================
#   Self Consistency test
#   (tests which neighborhood distance works best)
# =============================================================

    def self_consistency_test(self):
        if not self.function_proteins or self.graph is None:
            print("  ❌ No data for self-consistency test!")
            return 1 #can't test, use default n=1

        print("\n" + "=" * 60)
        print("📈 SELF-CONSISTENCY TEST")
        print("=" * 60)
        print("  🔍 Testing neighborhood distances (n=1,2,3):")
        print("-" * 60)

        # Get all proteins in the network
        network_proteins = list(self.graph.nodes())
        func_set = set(self.function_proteins)

        # Separate into functional and non-functional proteins
        functional_in_network = [p for p in network_proteins if p in func_set]
        non_functional_in_network = [p for p in network_proteins if p not in func_set]

        print("  📊 Found {} functional proteins in network".format(len(functional_in_network)))
        print("  📊 Found {} non-functional proteins in network".format(len(non_functional_in_network)))

        # Take equal numbers from both groups for testing
        test_size = min(25, len(functional_in_network), len(non_functional_in_network))
        if test_size < 5:
            print("  ⚠️ Not enough proteins for testing, using default n=1")
            return 1

        # create test set
        test_functional = functional_in_network[:test_size] #default: first 25 functional
        test_non_functional = non_functional_in_network[:test_size] #first 25 non-functional

        # Combine and create labels
        test_proteins = test_functional + test_non_functional
        is_functional = [True] * test_size + [False] * test_size
        print("  🧪 Testing on {} functional + {} non-functional proteins".format(test_size, test_size))

        # empty dictionary to store accuracies
        accuracies = {}
        # test each distance (n=1,2,3)
        for n in [1, 2, 3]:
            correct = 0
            tested = 0

            for i, protein in enumerate(test_proteins):
                if protein in self.graph:
                    # Get neighbors at distance n
                    neighbors = self.get_n_neighbors(protein, n)

                    # Count functional neighbors
                    annotated_neighbors = 0
                    for neighbor in neighbors:
                        if neighbor in func_set:  # Using func_set instead of self.function_proteins
                            annotated_neighbors += 1

                    if len(neighbors) > 0:
                        fraction = annotated_neighbors / len(neighbors)

                        # Predict; if >30% of neighbors are functional, it's functional
                        # (30% is more realistic than 50% for biological networks)
                        predict_has_function = fraction > 0.3
                        # Actual; does it really have the function?
                        actual_has_function = is_functional[i]

                        # check if prediction was correct
                        if predict_has_function == actual_has_function:
                            correct += 1
                        tested += 1

            if tested > 0:
                accuracy = correct / tested
                accuracies[n] = accuracy
                if accuracy >= 0.7:
                    print("  n={}: ✅ {:.1%}  ({}/{})".format(n, accuracy, correct, tested))
                elif accuracy >= 0.5:
                    print("  n={}: ⚠️ {:.1%}  ({}/{})".format(n, accuracy, correct, tested))
                else:
                    print("  n={}: ❌ {:.1%}  ({}/{})".format(n, accuracy, correct, tested))

        # Choose Best n
        if accuracies:
            self.best_n = max(accuracies, key=accuracies.get) # n with highest accuracy
            print("\n" + "=" * 60)
            print("✅ OPTIMAL DISTANCE: n={} (accuracy: {:.1%})".format(self.best_n, accuracies[self.best_n]))
            print("=" * 60)
        else:
            self.best_n = 1
            print("\n  ✅ USING DEFAULT: n={}".format(self.best_n))

        return self.best_n
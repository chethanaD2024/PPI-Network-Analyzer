"""
Name: A.D.D.C. Wijethunge
Index: s16562
PPI Network Analyzer - GUI
A graphical interface program to analyze PPI networks using Hishikagi algorithm for functional predictions

"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from ppi_analyzer import PPIAnalyzer

# Set seaborn style for all plots
sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1.2)

class PPIAnalyzerGUI:

# ============================================================
#   Initialize GUI application
# =============================================================

    def __init__(self):
        #create main window
        self.root = tk.Tk()
        self.root.title("PPI Network Neighborhood Analyzer")
        self.root.geometry("1000x700")

        #add GUI icon
        try:
            icon = tk.PhotoImage(file='PNetAnalyzer.png')
            self.root.iconphoto(True, icon)
        except:
            pass

        # initialize analyzer object (do the calculations)
        self.analyzer = PPIAnalyzer()
        # variable to store hishigaki scores after calculation(starts as 0)
        self.current_scores = None

        #create all GUI widgets
        self.create_widgets()


# ============================================================
#   Create and arrange widgets
# =============================================================

    def create_widgets(self):
        # -----Main container-----
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # -----Header-----
        title_label = tk.Label(
            main_frame,
            text="PPI NETWORK ANALYZER",
            font=("Arial", 18, "bold"),
            foreground="white",
            background="#2c3e50",
            padx=20,
            pady=10
        )
        title_label.pack(fill=tk.X , pady=(0, 10))

        # -----Separator line-----
        ttk.Separator(main_frame).pack(fill=tk.X, pady=(0, 20))

        # -----File selection Frame-----
        file_frame = ttk.LabelFrame(main_frame, text="Select Files", padding="10")
        file_frame.pack(fill=tk.X, pady=(0, 20))

        # PPI network file selection
        ttk.Label(file_frame, text="PPI Network File:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        # create text entry box
        self.ppi_entry = ttk.Entry(file_frame, width=50)
        # place the text entry box
        self.ppi_entry.grid(row=0, column=1, padx=5, pady=5)
        #creat "Browse" button
        ttk.Button(file_frame, text="Browse", command=self.browse_ppi).grid(row=0, column=2, padx=5)

        # Function protein file selection
        ttk.Label(file_frame, text="Function Protein File:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.func_entry = ttk.Entry(file_frame, width=50)
        self.func_entry.grid(row=1, column=1, padx=5, pady=5)
        ttk.Button(file_frame, text="Browse", command=self.browse_func).grid(row=1, column=2, padx=5)


        # -----Button panel-----
        #create frame to hold buttons
        button_frame = ttk.Frame(main_frame)
        #place the frame, stretches horizontally with 20px below
        button_frame.pack(fill=tk.X, pady=(0, 20))

        #main action buttons "Load Data"
        self.load_btn = ttk.Button(button_frame, text="Load Data", command=self.load_data)
        self.load_btn.pack(side=tk.LEFT, padx=5)

        #create run analysis button
        self.analyzer_btn = ttk.Button(button_frame, text="Run Analysis", command=self.run_analysis, state=tk.DISABLED)
        self.analyzer_btn.pack(side=tk.LEFT, padx=5)

        #create other buttons
        ttk.Button(button_frame, text="Show Plots", command=self.show_plots).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Save Results", command=self.save_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear", command=self.clear_all).pack(side=tk.LEFT, padx=5)


        # -----Notebook for tab control-----
        #create tabbed interface
        notebook = ttk.Notebook(main_frame)
        notebook.pack(fill=tk.BOTH, expand=True)

        # Console tab
        console_frame = ttk.Frame(notebook)
        notebook.add(console_frame, text="Console")

        #console text area for displaying output
        self.console_text = tk.Text(console_frame, height=20)
        # scroll bar to control text area's vertical view
        scrollbar = ttk.Scrollbar(console_frame, command=self.console_text.yview)
        # connects text area to the scroll bar
        self.console_text.config(yscrollcommand=scrollbar.set)
        #placing the text area
        self.console_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        #placing scroll bar
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)


        # Results tab
        results_frame = ttk.Frame(notebook)
        notebook.add(results_frame, text="Results")

        #treeview for displaying results in table format
        self.tree = ttk.Treeview(results_frame, columns=("Rank", "Protein", "Score"), show="headings")
        # set column headers
        self.tree.heading("Rank", text="Rank")
        self.tree.heading("Protein", text="Protein")
        self.tree.heading("Score", text="Hishigaki Score")
        # set column widths
        self.tree.column("Protein", width=200)
        self.tree.column("Score", width=150, anchor="e")
        # create and attaches a scroll bar to table
        tree_scroll = ttk.Scrollbar(results_frame, command=self.tree.yview)
        self.tree.config(yscrollcommand=tree_scroll.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        # Status bar
        self.status_bar = ttk.Label(main_frame, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(fill=tk.X, pady=(10, 0))

        # Redirect standard output to print on console tab
        self.redirect_stdout()


# ============================================================
#   Redirect print statement to console text widget
# =============================================================

    def redirect_stdout(self):
        class TextRedirector:
            def __init__(self, text_widget):
                self.text_widget = text_widget # inner class will capture print output

            def write(self, string): # when something is printed;
                self.text_widget.insert(tk.END, string) # insert text at the end of the console
                self.text_widget.see(tk.END) # scroll to bottom
                self.text_widget.update_idletasks() # updates display

            def flush(self): # requires for print redirection, does nothing
                pass

        # replace normal print destination with user TextRedirector
        sys.stdout = TextRedirector(self.console_text) # print() goes to console tab


# ============================================================
#   File Browsing
# =============================================================

    # method to call when the browse button is clicked
    # For PPI network file
    def browse_ppi(self):
        # open the file browse window, and return selected file path (or empty string if cancelled)
        filename = filedialog.askopenfilename(
            title="Select PPI Network File",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
        )
        if filename: # if file was selected;
            self.ppi_entry.delete(0, tk.END) # clear the entry box
            self.ppi_entry.insert(0, filename) # insert the file path

    # For function protein file
    def browse_func(self):
        filename = filedialog.askopenfilename(
            title="Select Function Protein File",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            self.func_entry.delete(0, tk.END)
            self.func_entry.insert(0, filename)


# ============================================================
#   Update the Status Bar with a message and color coding
# =============================================================

    def update_status(self, message):
        # change the text in the status bar
        self.status_bar.config(text=message)

        #color code based on message content
        msg_lower = message.lower()
        if "error" in msg_lower or "fail" in msg_lower:
            color = "red"
        elif "loading" in msg_lower or "running" in msg_lower:
            color = "blue"
        elif "ready" in msg_lower or "complete" in msg_lower:
            color = "green"
        else:
            color = "black"

        # applies the color and update the display
        self.status_bar.config(foreground=color)
        self.root.update_idletasks()


# ============================================================
#   Embed Matplotlib figure into a Tkinter frame
# =============================================================

    def embed_plot(self, fig, parent_frame):
        # create canvas to display figure
        canvas = FigureCanvasTkAgg(fig, master=parent_frame)
        # renders the figure
        canvas.draw()
        # get tkinter widget from canvas, packs it to fill the parent frame
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        return canvas # not used, but kept for reference




# ============================================================
#   Method for "Load Data" button clicks
# =============================================================

    def load_data(self):
        #load ppi network and function protein data
        ppi_file = self.ppi_entry.get() # get file paths from entry box
        func_file = self.func_entry.get()

        # Error message for not selecting PPI file
        if not ppi_file:
            messagebox.showerror("Error", "Please select a PPI network file")
            return

        self.update_status("Loading data...") # update status
        self.load_btn.config(state=tk.DISABLED) # user can't click the button twice
        self.console_text.delete(1.0, tk.END) # clear the console

        try:
            print("=" * 50)
            print("Loading PPI network...")

            # load ppi network using analyzer
            graph = self.analyzer.load_ppi_network(ppi_file)

            if graph:
                # Gets network statistics (number of proteins, interactions)
                proteins, interactions = self.analyzer.get_network_stats()

                # Save degree distribution plot
                print("\nSaving degree distribution plot...")
                #create 'result' directory if it doesn't exist
                results_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
                os.makedirs(results_dir, exist_ok=True)
                #save plot to results directory
                self.analyzer.plot_degree_distribution(save_path=os.path.join(results_dir, 'degree_distribution.png'))

                # Load function proteins if file provided
                if func_file and os.path.exists(func_file):
                    print("\nLoading function proteins...")
                    self.analyzer.load_function_proteins(func_file)

                # Enable 'Run Analysis' button
                self.analyzer_btn.config(state=tk.NORMAL)
                #update status with loaded stats
                self.update_status(f"Loaded: {proteins} proteins, {interactions} interactions")

                #test with sample protein
                # show degree of 1st protein in the network (just for demo)
                if self.analyzer.graph:
                    sample = list(self.analyzer.graph.nodes())[0]
                    self.analyzer.calculate_degree(sample)
            else:
                # Error message when the Graph loading fails
                messagebox.showerror("Error", "Failed to load PPI network")
                self.update_status("Load failed")

        # catches any errors, print them, show pops
        except Exception as e:
            print(f"Error: {str(e)}")
            messagebox.showerror("Error", f"Failed to load: \n{str(e)}")
            self.update_status("Error")
        # re-enables load data button (whether success or failure)
        finally:
            self.load_btn.config(state=tk.NORMAL)



# ============================================================
#   Method for "Run Analysis" button clicks
# =============================================================

    def run_analysis(self):
        # check network load
        if self.analyzer.graph is None:
            messagebox.showerror("Error", "Please load PPI network first")
            return

        #update staus
        self.update_status("Running analysis...")
        #disable button
        self.analyzer_btn.config(state=tk.DISABLED)

        try:
            print("\n" + "=" * 50)
            print("RUNNING ANALYSIS") # print header to console

            # Self-consistency test to determine optimal neighborhood distance
            print("\n1. Self-consistency test...")
            best_n = self.analyzer.self_consistency_test()
            print(f"Optimal neighborhood distance: n={best_n}")

            if self.analyzer.function_proteins: # only continues if seed proteins were loaded
                # count annotated neighbors for the 1st protein
                print("\n2. Counting annotated neighbors...")
                sample = list(self.analyzer.graph.nodes())[0]
                self.analyzer.count_annotated_neighbors(sample)

                # Save annotation distribution plot
                print("\nSaving annotation distribution plot...")
                #create results directory if it doesn't exist
                results_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
                os.makedirs(results_dir, exist_ok=True)
                self.analyzer.plot_annotation_distribution(
                    remove_zeros=True,
                    save_path=os.path.join(results_dir, 'annotation_distribution.png')
                )

                # Calculates Hishigaki scores and store them
                print("\n3. Calculating Hishigaki scores...")
                self.current_scores = self.analyzer.calculate_hishigaki_scores()

                
                # Saves Hishigaki distribution plot (Histogram)
                print("\nSaving Hishigaki distribution plot...")
                # Save to results directory (one level up from src)
                results_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
                os.makedirs(results_dir, exist_ok=True)
                self.analyzer.plot_hishigaki_distribution(
                    self.current_scores,
                    remove_zeros=True,
                    save_path=os.path.join(results_dir, 'hishigaki_distribution.png')
                )

                # Save Hishigaki Boxplot
                print("\nSaving Hishigaki boxplot...")
                self.analyzer.plot_hishigaki_boxplot(
                    self.current_scores,
                    remove_zeros=True,
                    save_path=os.path.join(results_dir, 'hishigaki_boxplot.png')
                )


                # Show results in table and updates status
                self.show_results()
                self.update_status(f"Analysis complete. Best n={best_n}")
            else:
                # message if no function is loaded
                print("\nNo function proteins loaded")
                self.update_status("Analysis complete (partial)")

        except Exception as e: # exception handling for errors
            print(f"Error: {str(e)}")
            messagebox.showerror("Error", f"Analysis failed:\n{str(e)}")
            self.update_status("Error")
        finally: # re-enables 'Run Analysis' button
            self.analyzer_btn.config(state=tk.NORMAL)



# ============================================================
#   Method for "Show Plots" button clicks
# =============================================================

    def show_plots(self):
        #check if data is  loaded, shows info if not
        if not self.analyzer or not self.analyzer.graph:
            messagebox.showinfo("No Data", "Please load data first!")
            return

        try:
            print("\n" + "=" * 50)
            print("Loading plots from analyzer...") # print header

            #create new window for plots (popup)
            plot_window = tk.Toplevel(self.root)
            plot_window.title("PPI Analysis Plots")
            plot_window.geometry("1200x900")

            # Create notebooks inside the new window
            plot_notebook = ttk.Notebook(plot_window)
            plot_notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)


            # Tab 1: Degree Distribution - USING ANALYZER'S METHOD
            degree_frame = ttk.Frame(plot_notebook)
            plot_notebook.add(degree_frame, text="Degree Distribution")

            # get degree plot from analyzer, embed it in the tab, then close the original figure (keeps a copy in the GUI)
            fig1 = self.analyzer.create_degree_plot_figure()
            if fig1:
                self.embed_plot(fig1, degree_frame)
                plt.close(fig1)
                print("✓ Degree distribution loaded")
            else:
                tk.Label(degree_frame, text="No degree data available").pack(pady=50)


            # Tab 2: Annotation Distribution - USING ANALYZER'S METHOD
            if self.analyzer.function_proteins:
                annot_frame = ttk.Frame(plot_notebook)
                plot_notebook.add(annot_frame, text="Annotation Distribution")

                #get annotation plot from analyzer
                fig2 = self.analyzer.create_annotation_plot_figure(remove_zeros=True)
                if fig2:
                    self.embed_plot(fig2, annot_frame)
                    plt.close(fig2)
                    print("✓ Annotation distribution loaded")


            # Tab 3: Hishigaki Scores - USING ANALYZER'S METHOD
            if self.current_scores:
                hish_frame = ttk.Frame(plot_notebook)
                plot_notebook.add(hish_frame, text="Hishigaki Scores")

                #get Hishigaki plot from analyzer (histogram + boxplot together)
                fig3 = self.analyzer.create_hishigaki_plot_figure(
                    self.current_scores,
                    remove_zeros=True
                )
                if fig3:
                    self.embed_plot(fig3, hish_frame)
                    plt.close(fig3)
                    print("✓ Hishigaki distribution loaded")


                # Tab 4: Top Candidates
                top_frame = ttk.Frame(plot_notebook)
                plot_notebook.add(top_frame, text="Top Candidates")

                # Get top 10 candidate scores  and extract protein names, scores
                top_10 = sorted(self.current_scores.items(), key=lambda x: x[1], reverse=True)[:10]
                proteins = [p[0] for p in top_10]
                scores = [p[1] for p in top_10]

                # create bar chart of top 10 candidates
                fig4, ax = plt.subplots(figsize=(14, 6))
                bars = ax.bar(range(len(proteins)), scores, color='skyblue', edgecolor='black')
                ax.set_xlabel('Protein', fontsize=12)
                ax.set_ylabel('Hishigaki Score', fontsize=12)
                ax.set_title('Top 10 Candidate Proteins', fontsize=14)
                ax.set_xticks(range(len(proteins)))
                ax.set_xticklabels(proteins, rotation=45)
                ax.grid(True, alpha=0.3)

                # Add score numbers on top of the each bar
                for bar, score in zip(bars, scores):
                    height = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width() / 2, height,
                           f'{score:.3f}', ha='center', va='bottom', fontsize=9)

                # Embed the bar chart
                plt.tight_layout()
                self.embed_plot(fig4, top_frame)
                plt.close(fig4)
                print("✓ Top candidates plot created")

            # print successful message
            print("✓ All plots displayed successfully")

        except Exception as e: # error handling
            print(f"Error creating plots: {str(e)}")
            import traceback
            traceback.print_exc()


# ============================================================
#   Method for Display scores in Results table
# =============================================================

    def show_results(self):
        # returns if no scores
        if not self.current_scores:
            return

        # clears existing results (all rows)
        self.tree.delete(*self.tree.get_children())

        # sort scores in descending order
        sorted_scores = sorted(self.current_scores.items(), key=lambda x: x[1], reverse=True)

        # display top 50 results in table ('enumerate' gives rank starting at 1)
        for rank, (protein, score) in enumerate(sorted_scores[:50], start=1): #only top 50
            self.tree.insert("", tk.END, values=(rank, protein, f"{score:.4f}"))

        # prints top 5 predictions to console
        print(f"\nTop 5 predictions:")
        for rank, (protein, score) in enumerate(sorted_scores[:5], start=1):
            print(f"  {rank}. {protein}: {score:.4f}")


# ============================================================
#   Method for "Save Results" button clicks
# =============================================================

    def save_results(self):
        # show warning of no results
        if not self.current_scores:
            messagebox.showwarning("No Data", "No analysis results to save!")
            return

        # default saving location is results directory (creates if it doesn't exist)
        default_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
        os.makedirs(default_dir, exist_ok=True)

        # ask useer for save location ('Save As' dialogue, default '.csv' format)
        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            initialfile="hishigaki_results.csv",
            initialdir=default_dir
        )

        if filename: # if user chose a file
            try:
                # calls analyzer's save method
                self.analyzer.save_results(self.current_scores, filename)
                #updates status
                self.update_status(f"Saved to {os.path.basename(filename)}")
            except Exception as e: # error handling
                messagebox.showerror("Save Error", f"Could not save: {str(e)}")



# ============================================================
#   Method for "Clear" button clicks
# =============================================================

    def clear_all(self):
        #ask user to confirm Y/N popup
        if messagebox.askyesno("Confirm", "Clear all data and results?"):
            # create a new empty analyzer and clears scores
            self.analyzer = PPIAnalyzer()
            self.current_scores = None

            #clears all entry boxes, console and table
            self.ppi_entry.delete(0, tk.END)
            self.func_entry.delete(0, tk.END)
            self.console_text.delete(1.0, tk.END)
            self.tree.delete(*self.tree.get_children())

            #disable 'Run Analysis' button
            self.analyzer_btn.config(state=tk.DISABLED)
            #updates status
            self.update_status("Cleared")
            #prints message
            print("All data cleared")


# ============================================================
#   Star GUI  Application
#   ( mainloop() runs forever, waiting for user clicks)
# =============================================================

    def run(self):
        self.root.mainloop()



# ============================================================
#   Main function to launch GUI application
# =============================================================

def main():
    app = PPIAnalyzerGUI()
    app.run()


if __name__ == "__main__":
    main()
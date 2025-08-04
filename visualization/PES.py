####################################################################
# Script: PES.py
# Author: Cristina Berga, https://github.com/CristinaBerga/CompChem-Tools/
# Creation Date: July 2025
####################################################################
#
# This Python script is a graphical user interface (GUI) tool
# for plotting energy profiles from data contained in Excel files.
# It allows users to visualize and customize reaction pathways with
# a high degree of control over the plot's appearance, making it
# ideal for presenting computational chemistry results.
#
####################################################################
#
# Key Features:
#
# 1. Interactive GUI: Provides a user-friendly interface built
#    with Tkinter to select input files and configure plot options.
#
# 2. Excel Integration: Reads energy data from multiple sheets
#    within a single Excel file, allowing the plotting of several
#    reaction pathways on the same graph.
#
# 3. Extensive Customization: Users can adjust numerous plot
#    parameters, including:
#    - Y-axis energy type (ΔE, ΔH, ΔG).
#    - Plot dimensions (A4 or Half width).
#    - Colors, line styles, and line widths for each pathway.
#    - Label and value positioning for each point.
#    - Plateau width and spacing.
#
# 4. High-Quality Output: Plots are generated using Matplotlib,
#    ensuring professional-quality graphics. The resulting figures
#    can be saved in PNG or SVG formats, with optimized settings
#    for publication.
#
####################################################################
#
# Usage:
#
# - Run the script to open the GUI.
#
# - Select an Excel file containing your energy profile data.
#   The script expects sheets with a column for labels and a
#   column for energy values (ΔE, ΔH, or ΔG).
#   See PES_Example.xlsx
#
# - Configure the plot settings as needed and click "Plot Energy
#   Profile" to generate the visualization.
#
####################################################################

import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox, colorchooser
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'

class EnergyProfilePlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("Energy Profile Plotter")

        self.sheet_configs = {}
        self.setup_ui()

    def setup_ui(self):
        self.file_path = filedialog.askopenfilename(title="Select Excel File", filetypes=[("Excel files", "*.xlsx *.xls")])
        if not self.file_path:
            return

        try:
            self.sheet_names = pd.ExcelFile(self.file_path).sheet_names
        except Exception as e:
            messagebox.showerror("Error", f"Could not read Excel file: {e}")
            return

        self.options_frame = tk.Frame(self.root)
        self.options_frame.pack(pady=10)

        tk.Label(self.options_frame, text="Select Sheets to Plot:").grid(row=0, column=0, sticky="w")
        self.sheet_vars = []
        for i, name in enumerate(self.sheet_names):
            var = tk.IntVar()
            chk = tk.Checkbutton(self.options_frame, text=name, variable=var)
            chk.grid(row=0, column=i+1, sticky="w")
            self.sheet_vars.append((var, name))

        self.energy_type = tk.StringVar(value="ΔE")
        tk.Label(self.options_frame, text="Y-axis (Energy Type):").grid(row=1, column=0, sticky="w")
        for i, val in enumerate(["ΔE", "ΔH", "ΔG"]):
            tk.Radiobutton(self.options_frame, text=val, variable=self.energy_type, value=val).grid(row=1, column=i+1, sticky="w")

        self.show_x = tk.BooleanVar(value=False)
        tk.Checkbutton(self.options_frame, text="Show X Axis (Reaction Coordinate)", variable=self.show_x).grid(row=2, column=0, columnspan=3, sticky="w")

        self.label_position = tk.StringVar(value="above")
        tk.Label(self.options_frame, text="Label Position:").grid(row=3, column=0, sticky="w")
        tk.Radiobutton(self.options_frame, text="Above", variable=self.label_position, value="above").grid(row=3, column=1, sticky="w")
        tk.Radiobutton(self.options_frame, text="Below", variable=self.label_position, value="below").grid(row=3, column=2, sticky="w")

        self.value_position = tk.StringVar(value="below")
        tk.Label(self.options_frame, text="Value Position:").grid(row=4, column=0, sticky="w")
        tk.Radiobutton(self.options_frame, text="Beside Label", variable=self.value_position, value="beside").grid(row=4, column=1, sticky="w")
        tk.Radiobutton(self.options_frame, text="Below Point", variable=self.value_position, value="below").grid(row=4, column=2, sticky="w")

        tk.Label(self.options_frame, text="Plateau Width:").grid(row=5, column=0, sticky="w")
        self.plateau_width = tk.DoubleVar(value=10.0)
        tk.Spinbox(self.options_frame, textvariable=self.plateau_width, from_=0.1, to=100.0, increment=0.5, width=5, format="%.1f").grid(row=5, column=1, sticky="w")

        tk.Label(self.options_frame, text="Plateau Spacing:").grid(row=6, column=0, sticky="w")
        self.plateau_spacing = tk.DoubleVar(value=30.0)
        tk.Spinbox(self.options_frame, textvariable=self.plateau_spacing, from_=0.1, to=100.0, increment=0.5, width=5, format="%.1f").grid(row=6, column=1, sticky="w")

        tk.Label(self.options_frame, text="Line Width:").grid(row=7, column=0, sticky="w")
        self.line_width = tk.DoubleVar(value=2.0)
        tk.Spinbox(self.options_frame, textvariable=self.line_width, from_=0.1, to=10.0, increment=0.1, width=5, format="%.1f").grid(row=7, column=1, sticky="w")

        tk.Label(self.options_frame, text="Dotted Line Width:").grid(row=8, column=0, sticky="w")
        self.dotted_width = tk.DoubleVar(value=1.0)
        tk.Spinbox(self.options_frame, textvariable=self.dotted_width, from_=0.1, to=10.0, increment=0.1, width=5, format="%.1f").grid(row=8, column=1, sticky="w")

        tk.Label(self.options_frame, text="Figure Size:").grid(row=9, column=0, sticky="w")
        self.figsize_option = tk.StringVar(value="A4")
        tk.Radiobutton(self.options_frame, text="A4 Width", variable=self.figsize_option, value="A4").grid(row=9, column=1, sticky="w")
        tk.Radiobutton(self.options_frame, text="Half Width", variable=self.figsize_option, value="Half").grid(row=9, column=2, sticky="w")

        self.show_legend = tk.BooleanVar(value=False)
        tk.Checkbutton(self.options_frame, text="Show Legend", variable=self.show_legend).grid(row=10, column=0, columnspan=2, sticky="w")
        
        tk.Label(self.options_frame, text="Decimal Places:").grid(row=11, column=0, sticky="w")
        self.decimal_places = tk.IntVar(value=1) 
        tk.Spinbox(self.options_frame, textvariable=self.decimal_places, from_=0, to=5, increment=1, width=5).grid(row=11, column=1, sticky="w")

        tk.Label(self.options_frame, text="Y-axis Tick Interval:").grid(row=12, column=0, sticky="w")
        self.y_tick_interval = tk.DoubleVar(value=2.0)
        tk.Spinbox(self.options_frame, textvariable=self.y_tick_interval, from_=0.1, to=20.0, increment=0.1, width=5, format="%.1f").grid(row=12, column=1, sticky="w")


        self.config_frame = tk.LabelFrame(self.root, text="Per Sheet Configuration")
        self.config_frame.pack(pady=10)

        self.sheet_config_vars = {}
        for idx, (_, name) in enumerate(self.sheet_vars):
            row = idx
            color_btn = tk.Button(self.config_frame, text=f"{name} Line Color", command=lambda n=name: self.pick_color(n, 'line'))
            color_btn.grid(row=row, column=0)
            text_color_btn = tk.Button(self.config_frame, text=f"{name} Text Color", command=lambda n=name: self.pick_color(n, 'text'))
            text_color_btn.grid(row=row, column=1)
            style_label = tk.Label(self.config_frame, text="Line Style:")
            style_label.grid(row=row, column=2)
            style_var = tk.StringVar(value='dotted')
            style_menu = tk.OptionMenu(self.config_frame, style_var, 'solid', 'dotted', 'dashed', 'dashdot')
            style_menu.grid(row=row, column=3)
            self.sheet_config_vars[name] = {'line': '#000000', 'text': '#000000', 'style': style_var}

        tk.Button(self.root, text="Plot Energy Profile", command=self.plot).pack(pady=10)

    def pick_color(self, name, key):
        color_code = colorchooser.askcolor(title=f"Choose {key} color for {name}")[1]
        if color_code:
            self.sheet_config_vars[name][key] = color_code

    def plot(self):
        selected_sheets = [name for var, name in self.sheet_vars if var.get()]
        if not selected_sheets:
            messagebox.showerror("Error", "No sheets selected")
            return

        fig_width = 6.69 if self.figsize_option.get() == "A4" else 3.35
        fig, ax = plt.subplots(figsize=(fig_width, 3))

        line_styles = {
            'solid': '-',
            'dotted': ':',
            'dashed': '--',
            'dashdot': '-.'
        }

        legend_handles = []
        decimal_places = self.decimal_places.get() 
        y_tick_interval = self.y_tick_interval.get() 

        for idx, sheet in enumerate(selected_sheets):
            df = pd.read_excel(self.file_path, sheet_name=sheet)
            labels = df.iloc[:, 0].tolist()
            values = [round(float(val), decimal_places) for val in df[self.energy_type.get()].tolist()]

            line_color = self.sheet_config_vars[sheet]['line']
            text_color = self.sheet_config_vars[sheet]['text']
            line_style_key = self.sheet_config_vars[sheet]['style'].get()
            line_style = line_styles.get(line_style_key, ':')

            label_positions = []
            start_x = 0
            for i, val in enumerate(values):
                x0 = start_x
                x1 = start_x + self.plateau_width.get()

                ax.plot([x0, x1], [val, val], color=line_color, linestyle='-', linewidth=self.line_width.get(), zorder=2)

                label_positions.append(((x0 + x1) / 2, val))
                start_x = x1 + self.plateau_spacing.get()

                if i < len(values) - 1:
                    next_val = values[i + 1]
                    x2 = start_x
                    ax.plot([x1, x2], [val, next_val], color=line_color, linestyle=line_style, linewidth=self.dotted_width.get(), zorder=1)

            for i, label in enumerate(labels):
                label_x, label_y = label_positions[i]
                value_str = f"({values[i]:.{decimal_places}f})"
                offset = 0.5 if self.label_position.get() == "above" else -1.2
                
                if self.value_position.get() == "beside":
                    full_label_text = f"{label} {value_str}"
                    ax.text(label_x, label_y + offset, full_label_text,
                            ha='center', va='bottom' if offset > 0 else 'top',
                            color=text_color)
                else:
                    ax.text(label_x, label_y + offset, label,
                            ha='center', va='bottom' if offset > 0 else 'top',
                            color=text_color)
                    ax.text(label_x, label_y - 1.0, value_str, ha='center', va='top', color=text_color)

            legend_handles.append(plt.Line2D([0], [0], color=line_color, linestyle=line_style, linewidth=self.line_width.get(), label=sheet))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if not self.show_x.get():
            ax.set_xticks([])
            ax.set_xlabel("")
            ax.spines['bottom'].set_visible(False)
        else:
            ax.set_xticks([])
            ax.set_xlabel("Reaction coordinate")

        ax.set_ylabel(f"{self.energy_type.get()} (kcal/mol)")
        ax.yaxis.set_major_locator(plt.MultipleLocator(y_tick_interval)) 

        if self.show_legend.get():
            ax.legend(handles=legend_handles, loc='upper right', frameon=True, edgecolor='black', shadow=True)

        plt.tight_layout()
        plt.show()

        save = messagebox.askyesno("Save", "Do you want to save the plot?")
        if save:
            filetypes = [("PNG", "*.png"), ("SVG", "*.svg")]
            save_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=filetypes)
            if save_path:
                if save_path.endswith(".png"):
                    fig.savefig(save_path, transparent=True, dpi=300)
                elif save_path.endswith(".svg"):
                    fig.savefig(save_path, format='svg', dpi=300, transparent=False, bbox_inches='tight', pad_inches=0.05, metadata={'Creator': 'EnergyProfilePlotter'})

if __name__ == '__main__':
    root = tk.Tk()
    app = EnergyProfilePlotter(root)
    root.mainloop()

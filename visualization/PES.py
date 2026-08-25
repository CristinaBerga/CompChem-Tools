####################################################################
# Script: PES.py
# Author: Cristina Berga, https://github.com/CristinaBerga/CompChem-Tools/
# Creation Date: July 2025
# Last Update: August 2026
####################################################################
#
# This Python script is a graphical user interface (GUI) tool
# for plotting energy profiles from data contained in Excel files.
# It allows users to visualize and customize reaction pathways with
# a high degree of control over the plot's appearance, making it
# ideal for presenting computational chemistry results. The graphs 
# can be exported in ChemDraw (CDXML), as well as in PNG and SVG formats.
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
from matplotlib import font_manager as fm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import math

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['path.simplify'] = False
mpl.rcParams['pdf.compression'] = 0

class CDXMLBuilder:
    
    FONT_ID = 3

    def __init__(self, font_name="Arial"):
        self._next_id = 100000
        self._body = []
        self.font_name = font_name
        
        self._custom_colors = []
        self._color_map = {
            '#000000': 3,
            '#000': 3,
            'BLACK': 3,
            '#FFFFFF': 2,
            '#FFF': 2,
            'WHITE': 2
        }

    def _id(self):
        self._next_id += 1
        return self._next_id

    def _get_color_index(self, hex_color):
        if not hex_color or not isinstance(hex_color, str):
            return 3
        
        hex_clean = hex_color.strip().upper()
        if not hex_clean.startswith('#') and len(hex_clean) in [3, 6]:
            hex_clean = '#' + hex_clean

        if hex_clean in self._color_map:
            return self._color_map[hex_clean]
        
        if len(hex_clean) == 7:
            r = int(hex_clean[1:3], 16) / 255.0
            g = int(hex_clean[3:5], 16) / 255.0
            b = int(hex_clean[5:7], 16) / 255.0
            self._custom_colors.append((r, g, b))
            idx = len(self._custom_colors) + 3
            self._color_map[hex_clean] = idx
            return idx
        return 3

    def add_line(self, x1, y1, x2, y2, width=0.6, linestyle='solid', color='#000000'):
        graphic_id = self._id()
        arrow_id = self._id()
        color_idx = self._get_color_index(color)
        
        ls_lower = str(linestyle).lower().strip()
        
        if ls_lower in ['dotted', 'dot']:
            line_attr = ' LineType="Dashed" HashSpacing="1.5" DashLength="1.5"'
        elif ls_lower in ['dashed', 'dash']:
            line_attr = ' LineType="Dashed" HashSpacing="3.0" DashLength="4.0"'
        elif ls_lower in ['dashdot', 'dash-dot', 'dash_dot']:
            line_attr = ' LineType="Dashed" HashSpacing="5.0" DashLength="3.0"'
        else:
            line_attr = ''

        self._body.append(f'''<graphic
 id="{graphic_id}"
 SupersededBy="{arrow_id}"
 BoundingBox="{min(x1,x2):.2f} {min(y1,y2):.2f} {max(x1,x2):.2f} {max(y1,y2):.2f}"
 GraphicType="Line"
 color="{color_idx}"{line_attr}
/>
<arrow
 id="{arrow_id}"
 BoundingBox="{min(x1,x2):.2f} {min(y1,y2):.2f} {max(x1,x2):.2f} {max(y1,y2):.2f}"
 FillType="None"
 ArrowheadType="None"
 color="{color_idx}"
 Head3D="{x2:.2f} {y2:.2f} 0"
 Tail3D="{x1:.2f} {y1:.2f} 0"
 LineWidth="{width}"{line_attr}
/>''')

    def add_text(self, text, x, y, size=10, bold=False, justification="Center", rotation_degrees=0, color='#000000'):
        text_str = str(text) if text is not None else ""
        if len(text_str) == 1:
            x -= 3
            
        color_idx = self._get_color_index(color)
        face = 1 if bold else 0
        rot = ''
        if rotation_degrees:
            rot = f'\n RotationAngle="{int(round(rotation_degrees * 65536))}"'
        text_id = self._id()

        self._body.append(f'''<t
 id="{text_id}"
 p="{x:.2f} {y:.2f}"
 BoundingBox="0 0 110 110"
 CaptionJustification="{justification}"
 Justification="{justification}"{rot}
 LineHeight="auto"
 font="{self.FONT_ID}"
 size="{size}"
 color="{color_idx}"
><s font="{self.FONT_ID}" size="{size}" color="{color_idx}" face="{face}">{text_str}</s></t>''')

    def render(self, name, page_width, page_height):
        colortable_entries = [
            '<color r="1.000" g="1.000" b="1.000"/>',
            '<color r="0.000" g="0.000" b="0.000"/>',
        ]
        for r, g, b in self._custom_colors:
            colortable_entries.append(f'<color r="{r:.3f}" g="{g:.3f}" b="{b:.3f}"/>')
        
        colortable_str = "\n".join(colortable_entries)
        page_id = self._id()
        body = "\n\n".join(self._body)
        
        return f'''<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
<CDXML
 CreationProgram="PES.py"
 Name="{name}"
 BoundingBox="0 0 {page_width:.2f} {page_height:.2f}"
 Font="{self.FONT_ID}"
 LabelFont="{self.FONT_ID}"
 LabelSize="10"
 CaptionFont="{self.FONT_ID}"
 CaptionSize="10"
 FormulaFont="{self.FONT_ID}"
 FormulaSize="10"
 LineWidth="0.60"
 BoldWidth="2.01"
 BondLength="14.40"
 color="3"
 bgcolor="2"
>
<colortable>
{colortable_str}
</colortable>
<fonttable>
<font id="{self.FONT_ID}" charset="iso-8859-1" name="{self.font_name}" family="Swiss"/>
</fonttable>
<page
 id="{page_id}"
 BoundingBox="0 0 {page_width:.2f} {page_height:.2f}"
 HeaderPosition="36"
 FooterPosition="36"
 PrintTrimMarks="yes"
 HeightPages="1"
 WidthPages="1"
>
{body}
</page>
</CDXML>
'''

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
        tk.Spinbox(self.options_frame, textvariable=self.line_width, from_=0.1, to=100.0, increment=0.1, width=5, format="%.1f").grid(row=7, column=1, sticky="w")

        tk.Label(self.options_frame, text="Dotted Line Width:").grid(row=8, column=0, sticky="w")
        self.dotted_width = tk.DoubleVar(value=1.0)
        tk.Spinbox(self.options_frame, textvariable=self.dotted_width, from_=0.1, to=100.0, increment=0.1, width=5, format="%.1f").grid(row=8, column=1, sticky="w")

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
        self.y_tick_interval = tk.DoubleVar(value=5.0)
        tk.Spinbox(self.options_frame, textvariable=self.y_tick_interval, from_=0.1, to=20.0, increment=0.1, width=5, format="%.1f").grid(row=12, column=1, sticky="w")

        tk.Label(self.options_frame, text="Font Family:").grid(row=13, column=0, sticky="w")
        system_fonts = sorted(list(set([f.name for f in fm.fontManager.ttflist])))
        if not system_fonts:
            system_fonts = ["Arial", "Times New Roman", "sans-serif", "serif"]
        default_font = "Arial" if "Arial" in system_fonts else system_fonts[0]
        self.font_family = tk.StringVar(value=default_font)
        self.font_menu = tk.OptionMenu(self.options_frame, self.font_family, *system_fonts)
        self.font_menu.grid(row=13, column=1, sticky="w")

        tk.Label(self.options_frame, text="Labels:").grid(row=14, column=0, sticky="w")
        self.text_display = tk.StringVar(value="Yes")
        self.display_menu = tk.OptionMenu(self.options_frame, self.text_display, "Yes", "Names Only", "Values Only", "No")
        self.display_menu.grid(row=14, column=1, sticky="w")

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

    def export_cdxml(self, save_path, selected_sheets):
        font_selected = self.font_family.get() if hasattr(self, 'font_family') else "Arial"
        b = CDXMLBuilder(font_name=font_selected)
        
        scale_x = 0.95
        pw = self.plateau_width.get() * scale_x * 2.2
        ps = self.plateau_spacing.get() * scale_x * 0.9
        
        axis_x = 55.0
        plateau_start_x = axis_x + 30.0

        all_values = []
        parsed_data = {}
        stage_order = []

        def get_base_stage(label):
            lbl = str(label).strip()
            for sep in ['-', '_', ' ']:
                if sep in lbl:
                    parts = lbl.split(sep)
                    if len(parts[1]) <= 2:
                        return parts[0]
            return lbl

        for sheet in selected_sheets:
            df = pd.read_excel(self.file_path, sheet_name=sheet)
            
            clean_labels = []
            clean_values = []
            
            for idx, row in df.iterrows():
                label_val = row.iloc[0]
                energy_val = row[self.energy_type.get()]
                
                if pd.notna(energy_val):
                    try:
                        val_float = float(energy_val)
                        if not math.isnan(val_float) and not math.isinf(val_float):
                            lbl_str = str(label_val).strip() if pd.notna(label_val) else f"P{len(clean_values)+1}"
                            clean_labels.append(lbl_str)
                            clean_values.append(round(val_float, self.decimal_places.get()))
                            
                            base_stage = get_base_stage(lbl_str)
                            if base_stage not in stage_order:
                                stage_order.append(base_stage)
                    except (ValueError, TypeError):
                        continue

            if clean_values:
                parsed_data[sheet] = list(zip(clean_labels, clean_values))
                all_values.extend(clean_values)

        if not all_values:
            return

        stage_x_coords = {}
        curr_x = plateau_start_x
        for stage in stage_order:
            stage_x_coords[stage] = curr_x
            curr_x += pw + ps

        max_x_reach = curr_x

        vmin, vmax = min(all_values), max(all_values)
        value_span = vmax - vmin or 1.0
        y_scale = 220.0 / value_span
        y_base = 280.0

        def y_of(v):
            return y_base - (v - vmin) * y_scale

        for sheet in selected_sheets:
            if sheet not in parsed_data:
                continue

            items = parsed_data[sheet]
            line_color = self.sheet_config_vars[sheet]['line']
            text_color = self.sheet_config_vars[sheet]['text']
            line_style = self.sheet_config_vars[sheet]['style'].get()

            drawn_points = []

            for lbl, val in items:
                base_stage = get_base_stage(lbl)

                if base_stage in stage_x_coords:
                    x0 = stage_x_coords[base_stage]
                else:
                    x0 = plateau_start_x + len(drawn_points) * (pw + ps)

                x1 = x0 + pw
                y = y_of(val)
                drawn_points.append((x1, y, x0))

                b.add_line(x0, y, x1, y, width=self.line_width.get(), linestyle='solid', color=line_color)

                display_mode = self.text_display.get()
                if display_mode != "No":
                    lx = (x0 + x1) / 2
                    val_str = f"({val:.{self.decimal_places.get()}f})"
                    is_above = (self.label_position.get() == "above")
                    offset_name = -10 if is_above else 14
                    offset_val = 12 if is_above else 24

                    if display_mode == "Yes":
                        if self.value_position.get() == "beside":
                            b.add_text(f"{lbl} {val_str}", lx, y + offset_name, bold=True, color=text_color)
                        else:
                            b.add_text(lbl, lx, y + offset_name, bold=True, color=text_color)
                            b.add_text(val_str, lx, y + offset_val, color=text_color)
                    elif display_mode == "Names Only":
                        b.add_text(lbl, lx, y + offset_name, bold=True, color=text_color)
                    elif display_mode == "Values Only":
                        b.add_text(val_str, lx, y + offset_name, color=text_color)

            for i in range(len(drawn_points) - 1):
                prev_x1, prev_y, _ = drawn_points[i]
                _, next_y, next_x0 = drawn_points[i + 1]
                b.add_line(prev_x1, prev_y, next_x0, next_y, width=self.dotted_width.get(), linestyle=line_style, color=line_color)

        step = self.y_tick_interval.get()
        lo = math.floor(vmin / step) * step
        hi = math.ceil(vmax / step) * step
        
        ticks = []
        t = lo
        while t <= hi + step * 0.5:
            ticks.append(round(t, 6))
            t += step

        axis_top_y = y_of(hi)
        axis_bottom_y = y_of(lo)
        b.add_line(axis_x, axis_top_y, axis_x, axis_bottom_y, width=1.0, linestyle='solid', color='#000000')
        
        for tick_val in ticks:
            ty = y_of(tick_val)
            b.add_line(axis_x - 5, ty, axis_x, ty, width=1.0, linestyle='solid', color='#000000')
            b.add_text(f"{tick_val:.1f}", axis_x - 10, ty + 3, justification="Right", color='#000000')

        axis_title_x = axis_x - 35
        axis_title_y = (axis_top_y + axis_bottom_y) / 2
        b.add_text(f"{self.energy_type.get()} (kcal/mol)", axis_title_x, axis_title_y, rotation_degrees=270, color='#000000')

        page_width = max_x_reach + 40.0
        page_height = 360.0
        
        with open(save_path, 'w', encoding='utf-8') as fh:
            fh.write(b.render("reaction_path.cdxml", page_width, page_height))

    def plot(self):
        selected_sheets = [name for var, name in self.sheet_vars if var.get()]
        if not selected_sheets:
            messagebox.showerror("Error", "No sheets selected")
            return

        plt.rcParams['font.family'] = self.font_family.get()

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

            display_mode = self.text_display.get()
            if display_mode != "No":
                for i, label in enumerate(labels):
                    label_x, label_y = label_positions[i]
                    value_str = f"({values[i]:.{decimal_places}f})"
                    offset = 0.5 if self.label_position.get() == "above" else -1.2
                    
                    if display_mode == "Yes":
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
                            
                    elif display_mode == "Names Only":
                        ax.text(label_x, label_y + offset, label,
                                ha='center', va='bottom' if offset > 0 else 'top',
                                color=text_color)
                                
                    elif display_mode == "Values Only":
                        ax.text(label_x, label_y + offset, value_str,
                                ha='center', va='bottom' if offset > 0 else 'top',
                                color=text_color)

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
            filetypes = [
                ("ChemDraw CDXML (*.cdxml)", "*.cdxml"),
                ("PNG Image (*.png)", "*.png"),
                ("SVG Vector (*.svg)", "*.svg"),
            ]
            save_path = filedialog.asksaveasfilename(defaultextension=".cdxml", filetypes=filetypes)
            if save_path:
                if save_path.endswith(".cdxml"):
                    self.export_cdxml(save_path, selected_sheets)
                elif save_path.endswith(".png"):
                    fig.savefig(save_path, transparent=True, dpi=300)
                elif save_path.endswith(".svg"):
                    fig.savefig(save_path, format='svg', dpi=300, transparent=False, bbox_inches='tight', pad_inches=0.05, metadata={'Creator': 'Cristina Berga'})

if __name__ == '__main__':
    root = tk.Tk()
    app = EnergyProfilePlotter(root)
    root.mainloop()

import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
import os
import sys

# Ensure we can import the package modules even if run as a script
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    __package__ = "ESRF_ID10_SURF"

from .XRR.XRR import XRR
from .GID.GID import GID
from .cli import parse_scans
from tkinter import ttk

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("ESRF ID10 SURF Processing")
        self.geometry("600x500")

        # Variables
        self.file_path = tk.StringVar()
        self.saving_path = tk.StringVar()
        self.mode = tk.StringVar(value="xrr")
        self.save_dat_var = tk.BooleanVar()
        self.save_orso_var = tk.BooleanVar()
        
        # GID Variables
        self.gid_scan_str = tk.StringVar()
        self.gid_PX0_str = tk.StringVar(value='50')
        self.gid_PPD_str = tk.StringVar(value='160')
        self.gid_mythen_gap_str = tk.StringVar(value='100')
        self.gid_alpha_i_str = tk.StringVar(value='mu')


        # XRR variables
        self.xrr_scan_str = tk.StringVar()
        self.z_scan_str = tk.StringVar()
        self.sample_size_str = tk.StringVar(value="1.0")
        self.z_scan_str = tk.StringVar()
        self.sample_size_str = tk.StringVar(value="1.0")
        self.beam_size_str = tk.StringVar(value="10.0")
        self.sample_size_str = tk.StringVar(value="1.0")
        self.beam_size_str = tk.StringVar(value="10.0")
        self.alpha_i_str = tk.StringVar(value="mu")
        self.px0_var = tk.StringVar(value="404")
        self.py0_var = tk.StringVar(value="165")

        # UI Layout
        self.create_widgets()

    def create_widgets(self):
        # File Selection
        frame_file = tk.Frame(self)
        frame_file.pack(pady=10, fill=tk.X, padx=10)
        
        lbl_file = tk.Label(frame_file, text="Source file:")
        lbl_file.pack(side=tk.LEFT)
        
        entry_file = tk.Entry(frame_file, textvariable=self.file_path, width=50)
        entry_file.pack(side=tk.LEFT, padx=5)
        
        btn_browse = tk.Button(frame_file, text="Browse", command=self.browse_file)
        btn_browse.pack(side=tk.LEFT)


        frame_saving_dir = tk.Frame(self)
        frame_saving_dir.pack(pady=10, fill=tk.X, padx=10)

        lbl_saving_dir = tk.Label(frame_saving_dir, text="Saving dir:")
        lbl_saving_dir.pack(side=tk.LEFT)

        entry_saving_dir = tk.Entry(frame_saving_dir, textvariable=self.saving_path, width=50)
        entry_saving_dir.pack(side=tk.LEFT, padx=5)

        btn_dir_browse = tk.Button(frame_saving_dir, text="Browse", command=self.browse_folder)
        btn_dir_browse.pack(side=tk.LEFT)

        # Mode Selection
        frame_mode = tk.Frame(self)
        frame_mode.pack(pady=10, fill=tk.X, padx=10)
        
        tk.Label(frame_mode, text="Mode:").pack(side=tk.LEFT)
        tk.Radiobutton(frame_mode, text="XRR", variable=self.mode, value="xrr", command=self.update_ui_state).pack(side=tk.LEFT)
        tk.Radiobutton(frame_mode, text="GID", variable=self.mode, value="gid", command=self.update_ui_state).pack(side=tk.LEFT)

        # Scans Container (Swappable frames)
        self.frame_scans_container = tk.Frame(self)
        self.frame_scans_container.pack(pady=5, fill=tk.X, padx=10)
        
        # XRR Scans Frame
        self.frame_xrr_inputs = tk.Frame(self.frame_scans_container)
        
        tk.Label(self.frame_xrr_inputs, text="XRR Scans:").grid(row=0, column=0, sticky='e', padx=5, pady=2)
        tk.Entry(self.frame_xrr_inputs, textvariable=self.xrr_scan_str, width=20).grid(row=0, column=1, sticky='w', padx=5, pady=2)
        
        tk.Label(self.frame_xrr_inputs, text="Z-Scan (Opt):").grid(row=0, column=2, sticky='e', padx=5, pady=2)
        tk.Entry(self.frame_xrr_inputs, textvariable=self.z_scan_str, width=10).grid(row=0, column=3, sticky='w', padx=5, pady=2)


        tk.Label(self.frame_xrr_inputs, text="Alpha_i:").grid(row=1, column=0, sticky='e', padx=5, pady=2)
        combo_alpha = ttk.Combobox(self.frame_xrr_inputs, textvariable=self.alpha_i_str, values=["mu", "chi"], width=8, state="readonly")
        combo_alpha.grid(row=1, column=1, sticky='w', padx=5, pady=2)


        tk.Label(self.frame_xrr_inputs, text="Sample Size (cm):").grid(row=2, column=0, sticky='e', padx=5, pady=2)
        tk.Entry(self.frame_xrr_inputs, textvariable=self.sample_size_str, width=10).grid(row=2, column=1, sticky='w', padx=5, pady=2)

        tk.Label(self.frame_xrr_inputs, text="Beam Size (um):").grid(row=2, column=2, sticky='e', padx=5, pady=2)
        tk.Entry(self.frame_xrr_inputs, textvariable=self.beam_size_str, width=10).grid(row=2, column=3, sticky='w', padx=5, pady=2)

        
        tk.Label(self.frame_xrr_inputs, text="PX0:").grid(row=3, column=0, sticky='e', padx=5, pady=2)
        tk.Entry(self.frame_xrr_inputs, textvariable=self.px0_var, width=10).grid(row=3, column=1, sticky='w', padx=5, pady=2)
        
        tk.Label(self.frame_xrr_inputs, text="PY0:").grid(row=3, column=1, sticky='e', padx=5, pady=2)
        tk.Entry(self.frame_xrr_inputs, textvariable=self.py0_var, width=10).grid(row=3, column=2, sticky='w', padx=5, pady=2)

        tk.Button(self.frame_xrr_inputs, text="Pick DB", command=self.open_picker).grid(row=3, column=3, padx=5, pady=2)

        # GID Scans Frame
        self.frame_gid_inputs = tk.Frame(self.frame_scans_container)

        tk.Label(self.frame_gid_inputs, text='GID Scan:').grid(row=0, column=0, sticky='e', padx=5, pady=2)
        tk.Entry(self.frame_gid_inputs, textvariable=self.gid_scan_str, width=10).grid(row=0, column=1, sticky='w', padx=5, pady=2)

        tk.Label(self.frame_gid_inputs, text='PX0:').grid(row=1, column=0, sticky='e', padx=5, pady=2)
        tk.Entry(self.frame_gid_inputs, textvariable=self.gid_PX0_str, width=10).grid(row=1, column=1, sticky='w', padx=5, pady=2)

        tk.Label(self.frame_gid_inputs, text='PPD:').grid(row=1, column=2, sticky='e', padx=5, pady=2)
        tk.Entry(self.frame_gid_inputs, textvariable=self.gid_PPD_str, width=10).grid(row=1, column=3, sticky='w', padx=5, pady=2)

        tk.Label(self.frame_gid_inputs, text='Gap in pixels:').grid(row=2, column=0, sticky='e', padx=5, pady=2)
        tk.Entry(self.frame_gid_inputs, textvariable=self.gid_mythen_gap_str, width=10).grid(row=2, column=1, sticky='w', padx=5, pady=2)

        
        # Options
        frame_options = tk.Frame(self)
        frame_options.pack(pady=10, fill=tk.X, padx=10)
        
        self.chk_save_dat = tk.Checkbutton(frame_options, text="Save .dat/.h5", variable=self.save_dat_var)
        self.chk_save_dat.pack(side=tk.LEFT)
        
        self.chk_save_orso = tk.Checkbutton(frame_options, text="Save ORSO", variable=self.save_orso_var)
        self.chk_save_orso.pack(side=tk.LEFT)

        # Initialize UI State
        self.update_ui_state()

        # Process Button
        tk.Button(self, text="Process", command=self.process, bg="lightblue", height=2).pack(pady=20, fill=tk.X, padx=50)

    def browse_file(self):
        filename = filedialog.askopenfilename(filetypes=[("HDF5 files", "*.h5"), ("All files", "*.*")], initialdir="/", title="Select .h5 File")
        if filename:
            self.file_path.set(filename)

    def browse_folder(self):
        folder = filedialog.askdirectory(initialdir="/", title="Select saving directory")
        if folder:
            self.saving_path.set(folder)

    def update_ui_state(self):
        mode = self.mode.get()
        
        # Toggle Scan Frames
        self.frame_xrr_inputs.pack_forget()
        self.frame_gid_inputs.pack_forget()
        
        if mode == "xrr":
            self.frame_xrr_inputs.pack(fill=tk.X)
            self.chk_save_orso.config(state="normal")
        else:
            self.frame_gid_inputs.pack(fill=tk.X)
            self.chk_save_orso.config(state="disabled")
            self.save_orso_var.set(False)

    def process(self):
        filename = self.file_path.get()
        saving_path = self.saving_path.get()
        if not filename:
            messagebox.showerror("Error", "Please select a file.")
            return

        if not saving_path:
            messagebox.showerror("Error", "Please select a saving directory.")
            return
        
        mode = self.mode.get()
        
        if mode == 'xrr':
             self.process_xrr_flow(filename, saving_path)
        elif mode == 'gid':
             self.process_gid_flow(filename, saving_path)

    def process_xrr_flow(self, filename, saving_path):
        scan_input = self.xrr_scan_str.get()
        if not scan_input:
            messagebox.showerror("Error", "Please enter XRR scan numbers.")
            return
            
        try:
            scans = parse_scans(scan_input)
        except ValueError:
            messagebox.showerror("Error", "Invalid XRR scan format.")
            return
            
        # Z-Scan (Optional)
        z_scan_input = self.z_scan_str.get()
        z_processor = None
        if z_scan_input:
            try:
                z_scans = parse_scans(z_scan_input)
                z_processor = XRR(filename, z_scans, alpha_i_name='zgH')
            except Exception as e:
                messagebox.showwarning("Warning", f"Could not load Z-scan: {e}\nProcessing without normalization.")
        
        # Dimensions
        try:
            sample_size = float(self.sample_size_str.get())
            beam_size = float(self.beam_size_str.get())
        except ValueError:
             messagebox.showerror("Error", "Sample size and Beam size must be numbers.")
             return

        # ORSO Metadata
        orso_metadata = {}
        if self.save_orso_var.get():
            orso_metadata = self.ask_orso_metadata()
            if orso_metadata is None: return

        try:
            alpha_i = self.alpha_i_str.get()
            if alpha_i == "mu":
                 monitor = "ionch2"
                 
            else:
                 monitor = "mon"
            px0 = int(self.px0_var.get())
            py0 = int(self.py0_var.get())
        except ValueError:
             messagebox.showerror("Error", "PX0 and PY0 must be integers.")
             return

        # Z-Scan (Optional)
        z_scan_input = self.z_scan_str.get()
        z_processor = None
        if z_scan_input:
            try:
                z_scans = parse_scans(z_scan_input)
                z_processor = XRR(filename, z_scans, alpha_i_name='zgH', monitor_name=monitor, PX0=px0, PY0=py0)
            except Exception as e:
                messagebox.showwarning("Warning", f"Could not load Z-scan: {e}\nProcessing without normalization.")
       
        try:
            processor = XRR(filename, scans, alpha_i_name=alpha_i, monitor_name=monitor, PX0=px0, PY0=py0, saving_dir=saving_path)
            
            # Auto corrections if params provided
            if z_processor:
                 processor.apply_auto_corrections(sample_size=sample_size, beam_size=beam_size, z_scan=z_processor)
            
            if self.save_dat_var.get():
                processor.save_reflectivity(format='dat')
            
            if self.save_orso_var.get():
                processor.save_reflectivity(format='orso', owner=orso_metadata['owner'], creator=orso_metadata['creator'])

            self.show_graph(processor, "xrr")
            
        except Exception as e:
            messagebox.showerror("Processing Error", str(e))

    def process_gid_flow(self, filename, saving_path):
        scan_input = self.gid_scan_str.get()
        if not scan_input:
            messagebox.showerror("Error", "Please enter GID scan numbers.")
            return

        try:
            scans = parse_scans(scan_input)
        except ValueError:
            messagebox.showerror("Error", "Invalid GID scan format.")
            return

        try:
            gid_px0 = float(self.gid_PX0_str.get())
            gid_ppd = float(self.gid_PPD_str.get())
            gid_mythen_gap = int(self.gid_mythen_gap_str.get())
        except ValueError:
            messagebox.showerror("Error", "Fill all experimental parameters.")
            return


        try:
            processor = GID(filename, scans, PX0=gid_px0, PPD=gid_ppd, mythen_gap=gid_mythen_gap, save_path=saving_path)
            
            if self.save_dat_var.get():
                processor.save_image_h5()
            
            self.show_graph(processor, "gid")
        except Exception as e:
            messagebox.showerror("Processing Error", str(e))



    def open_picker(self):
        filename = self.file_path.get()
        scan_input = self.xrr_scan_str.get()
        if not filename or not scan_input:
             messagebox.showerror("Error", "Please select file and XRR scans first.")
             return
             
        try:
            scans = parse_scans(scan_input)
            if not scans: return
            first_scan = [scans[0]]

            px0 = int(self.px0_var.get())
            py0 = int(self.py0_var.get())
            
            # Initialize minimal XRR to get data
            tmp_processor = XRR(filename, first_scan, alpha_i_name=self.alpha_i_str.get(), PX0=px0, PY0=py0)
            
            top = tk.Toplevel(self)
            top.title(f"Pick Direct Beam (Scan {first_scan[0]})")
            
            fig = plt.Figure(figsize=(6, 6), dpi=100)
            ax = fig.add_subplot(111)
            
            try:
                n_frames = len(tmp_processor.data)
                frame = min(5, n_frames // 2)
                tmp_processor.show_detector_image(frame_number=frame, ax=ax, plot_cross=True)
                ax.set_xlim(0, max(tmp_processor.data.shape))
                ax.set_ylim(0, max(tmp_processor.data.shape))
                ax.set_title(f"Click on Direct Beam Center (Scan {first_scan[0]}, Frame {frame})")
            except Exception as e:
                ax.text(0.5, 0.5, f"Error showing image: {e}", ha='center')
            
            canvas = FigureCanvasTkAgg(fig, master=top)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

            toolbar = NavigationToolbar2Tk(canvas, top)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            
            def onclick(event):
                if toolbar.mode != "": return
                if event.inaxes == ax:
                    ix, iy = event.xdata, event.ydata
                    if ix is not None and iy is not None:
                        self.px0_var.set(str(int(ix)))
                        self.py0_var.set(str(int(iy)))
                        messagebox.showinfo("Selected", f"Selected (PX0, PY0) = ({int(ix)}, {int(iy)})")
                        top.destroy()
            
            cid = fig.canvas.mpl_connect('button_press_event', onclick)
            
        except Exception as e:
            messagebox.showerror("Picker Error", str(e))

    def ask_orso_metadata(self):
        dialog = OrsoDialog(self)
        self.wait_window(dialog.top)
        return dialog.result



    def show_graph(self, processor, mode):
        top = tk.Toplevel(self)
        top.title(f"{mode.upper()} Result")
        top.geometry("800x800")

        fig = plt.Figure(figsize=(6, 6), dpi=100)

        if mode == 'xrr':
            ax = fig.add_subplot(111)
            processor.plot_reflectivity(ax=ax)
        elif mode == 'gid':
             processor.plot_quick_analysis(fig=fig)


        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, top)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        btn_save_graph = tk.Button(top, text="Save Graph", command=lambda: self.save_graph(fig))
        btn_save_graph.pack(pady=10)

        plt.close()

    def save_graph(self, fig):
        filename = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png"), ("PDF files", "*.pdf")],initialdir=self.saving_path.get())
        if filename:
            fig.savefig(filename)
            messagebox.showinfo("Saved", f"Graph saved to {filename}")


class OrsoDialog:
    def __init__(self, parent):
        self.top = tk.Toplevel(parent)
        self.top.title("ORSO Metadata")
        self.result = None

        tk.Label(self.top, text="Owner:").grid(row=0, column=0, padx=5, pady=5)
        self.owner_var = tk.StringVar(value="User University")
        tk.Entry(self.top, textvariable=self.owner_var).grid(row=0, column=1, padx=5, pady=5)

        tk.Label(self.top, text="Creator:").grid(row=1, column=0, padx=5, pady=5)
        self.creator_var = tk.StringVar(value="User Name")
        tk.Entry(self.top, textvariable=self.creator_var).grid(row=1, column=1, padx=5, pady=5)

        tk.Button(self.top, text="OK", command=self.on_ok).grid(row=2, column=0, columnspan=2, pady=10)

    def on_ok(self):
        self.result = {
            "owner": self.owner_var.get(),
            "creator": self.creator_var.get()
        }
        self.top.destroy()

if __name__ == "__main__":
    app = App()
    app.mainloop()

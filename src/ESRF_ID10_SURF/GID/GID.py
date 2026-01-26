import h5py
import time
import os
import copy
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import lmfit
from lmfit.models import GaussianModel, LorentzianModel, VoigtModel, PseudoVoigtModel, LinearModel, ConstantModel

from ..base import BaseSurf

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

# --- Main GID Class ---

class GID(BaseSurf):
    """
    Main class for processing GID (Grazing Incidence Diffraction) data.
    """

    def __init__(self, file, scans, alpha_i_name='chi', detector_name='mythen2', monitor_name='mon',
                 transmission_name='autof_eh1_transm', att_name='autof_eh1_curratt', cnttime_name='sec',
                 PX0=50, mythen_gap=120, PPD=198.5, pixel_size_qxz=0.055, angle_name='delta', energy_name='monoe',
                 I0=1e12, saving_dir=None, *args, **kwargs):

        super().__init__(file, scans, alpha_i_name=alpha_i_name, detector_name=detector_name,
                         monitor_name=monitor_name, transmission_name=transmission_name,
                         att_name=att_name, cnttime_name=cnttime_name, energy_name=energy_name,
                         I0=I0, saving_dir=saving_dir, **kwargs)

        self.angle_name = angle_name
        self.PX0 = PX0
        self.mythen_gap = mythen_gap
        self.PPD = PPD
        self.pixel_size_qxz = pixel_size_qxz

        # Attributes specific to GID
        self.angle = np.array([])

        # Processed data containers
        self.data_gap = None
        self.data_gap_e = None
        self.qz = None
        self.qxy = None

        # Add GID specific appendable keys
        self.appendable_keys.extend(['data', 'angle'])

        self.load_data()
        self.__process_2D_data__()
        self._check_saving_dir()

    def _load_single_scan(self, ScanN):
        """
        Load data for a single scan from the HDF5 file.
        """
        logger.info('Loading scan #{}'.format(ScanN))
        data_dict = {}
        try:
            with h5py.File(self.file, "r") as f:
                # Using [()] to read dataset into numpy array immediately
                data_dict['data'] = f.get(f"{ScanN}.1/measurement/{self.detector_name}")[()]
                data_dict['angle'] = f.get(f"{ScanN}.1/measurement/{self.angle_name}")[()]
                data_dict['alpha_i'] = f.get(f"{ScanN}.1/instrument/positioners/{self.alpha_i_name}")[()]
                data_dict['monitor'] = f.get(f"{ScanN}.1/measurement/{self.monitor_name}")[()]
                data_dict['transmission'] = f.get(f"{ScanN}.1/measurement/{self.transmission_name}")[()]
                data_dict['attenuator'] = f.get(f"{ScanN}.1/measurement/{self.att_name}")[()]
                data_dict['cnttime'] = f.get(f"{ScanN}.1/measurement/{self.cnttime_name}")[()]

                energy = f.get(f"{ScanN}.1/instrument/positioners/{self.energy_name}")[()]
                data_dict['energy'] = float(energy)

                sample_name_ds = f.get(f"{ScanN}.1/sample/name/")
                if sample_name_ds:
                    # Handle string decoding if necessary
                    data_dict['sample_name'] = str(sample_name_ds[()])[2:-1:1]

                pi_ds = f.get(f"{ScanN}.1/measurement/fb_Pi")
                if pi_ds:
                    Pi = np.mean(pi_ds[()])
                    if Pi < 90:
                        data_dict['Pi'] = int(np.round(Pi, 0))
                    else:
                         data_dict['Pi'] = '' # Or some default
        except Exception as e:
            logger.error(f"Error loading scan {ScanN}: {e}")
            raise

        return data_dict

    def get_qz(self, pixels):
        """
        Calculate the vertical scattering vector qz.
        """
        # Calculate qz. Assuming alpha_i is scalar or matches dimensions if it varies.
        # Here we treat alpha_i as a scalar (first value) if it's an array to produce a 1D qz array for pixels.
        alpha_i = self.alpha_i[0] if np.size(self.alpha_i) > 1 else self.alpha_i

        # Ensure energy is scalar
        energy = self.energy if np.ndim(self.energy) == 0 else self.energy[0]

        wavelength = 12.398 / energy   # in Angstroms
        k0 = 2 * np.pi / wavelength

        # pixels are an array
        qz = k0 * (np.sin(np.deg2rad(alpha_i)) + np.sin(np.deg2rad((pixels - self.PX0) / self.PPD)))
        return qz

    def get_qxy(self, angle):
        """
        Calculate the horizontal scattering vector qxy.
        """
        energy = self.energy if np.ndim(self.energy) == 0 else self.energy[0]
        wavelength = 12.398 / energy
        k0 = 2 * np.pi / wavelength
        qxy = 2 * k0 * np.sin(np.deg2rad(angle / 2))
        return qxy

    def __process_2D_data__(self):
        """
        Process the raw 2D detector data.
        """
        t0 = time.time()
        logger.info("Start processing 2D data.")
        nx, ny = np.shape(self.data)

        # Handle gaps in detector modules (Mythen is built from 2 modules)
        map2Dm = np.ones((nx, ny + self.mythen_gap))

        if ny >= 2559:
            map2Dm[:, 0:1279] = self.data[:, 0:1279]
            map2Dm[:, (1280 + self.mythen_gap):(2559 + self.mythen_gap)] = self.data[:, 1280:2559]
        else:
            logger.warning(f"Warning: Data shape {self.data.shape} does not match expected Mythen format. Using raw data.")
            map2Dm = self.data

        nxm, nym = np.shape(map2Dm)

        norm_factor = self.transmission * self.monitor / self.monitor[0]
        # Avoid division by zero
        norm_factor = np.where(norm_factor == 0, 1.0, norm_factor)

        self.data_gap = map2Dm / norm_factor[:, np.newaxis] / self.I0
        self.data_gap_e = np.sqrt(map2Dm) / norm_factor[:, np.newaxis] / self.I0

        self.qz = self.get_qz(np.arange(nym))
        self.qxy = self.get_qxy(self.angle)

        logger.info("Processing completed. Processing time %3.3f sec \n\n" % (time.time() - t0))

    def plot_2D_image(self, ax=None, save=False, **kwargs):
        """
        Plot the 2D GID map in reciprocal space (qxy vs qz).
        """
        if ax is None:
            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), layout='tight')
        else:
            ax0 = ax

        mean_val = np.mean(self.data_gap)
        std_val = np.std(self.data_gap)
        _vmin = np.log10(max(mean_val - std_val, 1e-12))
        _vmax = np.log10(mean_val + 3 * std_val)

        im = ax0.imshow(np.log10(np.rot90(self.data_gap)), aspect='equal', vmin=_vmin, vmax=_vmax,
                        extent=(np.min(self.qxy), np.max(self.qxy), np.min(self.qz), np.max(self.qz)), **kwargs)

        ax0.set_xlabel('$q_{xy}, \\AA^{-1}$')
        ax0.set_ylabel('$q_{z}, \\AA^{-1}$')

        if save:
            logger.info('Saving 2D GID map.')
            self._save_figure(plt.gcf(), '2Dmap')

        return ax0

    def get_qxy_cut(self, qz_min, qz_max):
        """
        Extract a 1D cut along qxy by integrating over a range of qz.
        """
        qz_indices = np.where((self.qz > qz_min) & (self.qz < qz_max))[0]
        if len(qz_indices) == 0:
            logger.warning(f"Warning: No data in qz range [{qz_min}, {qz_max}]")
            return [self.qxy, np.zeros_like(self.qxy)]

        cut_qxy = np.sum(self.data_gap[:, qz_indices], axis=1)
        return [self.qxy, cut_qxy]

    def save_qxy_cut(self, qz_min, qz_max, **kwargs):
        """
        Save the qxy cut data to a text file.
        """
        qxy, cut_qxy = self.get_qxy_cut(qz_min, qz_max)
        out = np.array([qxy, cut_qxy])

        self._ensure_sample_dir()

        filename = self.saving_dir + '/GID_{}_scan_{}_qxy_cut_{}_{}_A.dat'.format(
            self.sample_name, self.scans, qz_min, qz_max)

        np.savetxt(filename, out.T)
        logger.info('GID cut saved as: {}'.format(filename))

    def _plot_cut(self, x, y, xlabel, ylabel, label, ax=None, save=False, filename=None, **kwargs):
        """
        Generic plotting function for 1D cuts.
        """
        if ax is None:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), layout='tight')

        ax.plot(x, y, 'o', markersize=5, label=label, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()

        if save:
            if filename is None:
                raise ValueError("Filename must be provided if save is True.")
            logger.info(f'Saving plot to {filename}')
            plt.savefig(filename)

    def plot_qxy_cut(self, qz_min, qz_max, ax=None, save=False, **kwargs):
        """
        Plot the qxy cut.
        """
        qxy, cut_qxy = self.get_qxy_cut(qz_min, qz_max)
        label = f'$Cut\\: {qz_min:.2f}<q_z<{qz_max:.2f}$'
        filename = self.saving_dir + f'/qxy_cut_{qz_min}_{qz_max}_A.png' if save else None
        self._plot_cut(qxy, cut_qxy, '$q_{xy}, \\AA^{-1}$', 'Intensity', label, ax, save, filename, **kwargs)

    def get_qz_cut(self, qxy_min, qxy_max):
        """
        Extract a 1D cut along qz by integrating over a range of qxy.
        """
        qxy_indices = np.where((self.qxy > qxy_min) & (self.qxy < qxy_max))[0]
        if len(qxy_indices) == 0:
            logger.warning(f"Warning: No data in qxy range [{qxy_min}, {qxy_max}]")
            return [self.qz, np.zeros_like(self.qz)]

        cut_qz = np.sum(self.data_gap[qxy_indices, :], axis=0)
        return [self.qz, cut_qz]

    def save_qz_cut(self, qxy_min, qxy_max, **kwargs):
        """
        Save the qz cut data to a text file.
        """
        qz, cut_qz = self.get_qz_cut(qxy_min, qxy_max)
        out = np.array([qz, cut_qz])

        self._ensure_sample_dir()

        filename = self.saving_dir + '/GID_{}_scan_{}_qz_cut_{}_{}_A.dat'.format(
            self.sample_name, self.scans, qxy_min, qxy_max)

        np.savetxt(filename, out.T)
        logger.info('GID cut saved as: {}'.format(filename))

    def plot_qz_cut(self, qxy_min, qxy_max, ax=None, save=False, **kwargs):
        """
        Plot the qz cut.
        """
        qz, cut_qz = self.get_qz_cut(qxy_min, qxy_max)
        label = f'$Cut\\: {qxy_min:.2f}<q_{{xy}}<{qxy_max:.2f}$'
        filename = self.saving_dir + f'/qz_cut_{qxy_min}_{qxy_max}_A.png' if save else None
        self._plot_cut(qz, cut_qz, '$q_{z}, \\AA^{-1}$', 'Intensity', label, ax, save, filename, **kwargs)

    def plot_quick_analysis(self, save=False, fig=None):
        """
        Perform a quick standard analysis plot including the 2D map and representative qxy cuts.
        """
        if fig is None:
            fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6, 6), layout='tight')
        else:
            ax = fig.subplots(1,2)

        self.plot_2D_image(ax=ax[0])

        # Default limits, might need adjustment based on data
        ax[0].set_ylim(0, 1.5)
        ax[0].set_xlim(np.min(self.qxy), np.max(self.qxy))

        # Cut lines
        ax[0].hlines([0.01, 0.3], 1.2, 1.6, linestyle='--', alpha=0.8, color='C0')
        ax[0].hlines([0.5, 0.95], 1.2, 1.6, linestyle='--', alpha=0.8, color='C1')

        self.plot_qxy_cut(0.01, 0.3, ax=ax[1])
        self.plot_qxy_cut(0.5, 0.95, ax=ax[1])

        plt.suptitle('GID : Sample {}, Scan {}, Pi = {} mN/m'.format(self.sample_name, self.scans, self.Pi))

        if save:
            logger.info('Saving standard GID plot.')
            self._save_figure(fig, 'quick_analysis')

        return fig, ax

    # --- New Features ---

    def fit_profile(self, x, y, model='gaussian', background='linear', limits=None, **kwargs):
        """
        Fit a profile to the specified model with background using lmfit.
        """

        # Apply limits
        if limits:
            mask = (x >= limits[0]) & (x <= limits[1])
            x_fit = x[mask]
            y_fit = y[mask]
        else:
            x_fit = x
            y_fit = y

        # Select model
        if model == 'gaussian':
            peak = GaussianModel(prefix='peak_')
        elif model == 'lorentzian':
            peak = LorentzianModel(prefix='peak_')
        elif model == 'voigt':
            peak = VoigtModel(prefix='peak_')
        elif model == 'pseudo_voigt':
            peak = PseudoVoigtModel(prefix='peak_')
        else:
            raise ValueError(f"Unknown model: {model}")

        # Select background
        if background == 'linear':
            bg = LinearModel(prefix='bg_')
        elif background == 'constant':
            bg = ConstantModel(prefix='bg_')
        else:
            bg = None

        if bg:
            mod = peak + bg
        else:
            mod = peak

        params = mod.make_params()

        # Initial guesses
        if background == 'linear':
            slope = (y_fit[-1] - y_fit[0]) / (x_fit[-1] - x_fit[0])
            intercept = y_fit[0] - slope * x_fit[0]
            params['bg_slope'].set(value=slope)
            params['bg_intercept'].set(value=intercept)
        elif background == 'constant':
            params['bg_c'].set(value=np.min(y_fit))

        peak_params = peak.guess(y_fit, x=x_fit)
        params.update(peak_params)

        result = mod.fit(y_fit, params, x=x_fit, **kwargs)

        return result

    def analyze_peak(self, x, y, model='voigt', background='linear', limits=None, save=False, filename_prefix='fit_result', **kwargs):
        """
        Wrapper for fit_profile to perform analysis, plotting, and saving.
        """
        logger.info(f"Fitting {model} profile...")
        result = self.fit_profile(x, y, model=model, background=background, limits=limits, **kwargs)

        # Plot
        fig, ax = plt.subplots(1,1,figsize=(6,6), layout='tight')

        ax.plot(x, y, 'o', label='Data', markersize=4, alpha=0.6)

        if limits:
            mask = (x >= limits[0]) & (x <= limits[1])
            x_fit = x[mask]
            ax.plot(x_fit, result.best_fit, 'r-', label='Fit', linewidth=3)
            ax.axvline(limits[0], color='k', linestyle='--', alpha=0.3)
            ax.axvline(limits[1], color='k', linestyle='--', alpha=0.3)
            fit_line = np.array([x_fit, result.best_fit])
        else:
            ax.plot(x, result.best_fit, 'r-', label='Fit', linewidth=2)
            fit_line = np.array([x, result.best_fit])

        ax.legend()
        ax.set_xlabel('q')
        ax.set_ylabel('Intensity')
        ax.set_title(f'{model.capitalize()} Peak Fit')

        if save:
            self._ensure_sample_dir()
            fname_base = f"{self.saving_dir}/{filename_prefix}_{self.sample_name}"

            fig_name = f"{fname_base}.png"
            fig.savefig(fig_name, dpi=100)
            logger.info(f"Graph saved to {fig_name}")

            txt_name = f"{fname_base}.txt"
            with open(txt_name, 'w') as f:
                f.write(result.fit_report())
                f.write('\n __________________________\n')
                np.savetxt(f, fit_line.T)
            logger.info(f"Fit parameters saved to {txt_name}")

        return result

    def save_image_h5(self, filename=None):
        """
        Save the 2D image data to an HDF5 file in q-coordinates.
        """
        if filename is None:
            self._ensure_sample_dir()
            filename = self.saving_dir +'/GID_{}_2D.h5'.format(
                self.sample_name, self.scans)

        try:
            with h5py.File(filename, 'a') as hf:
                scan = hf.create_group(str(self.scans))

                scan.create_dataset('intensity', data=self.data_gap)
                scan.create_dataset('qxy', data=self.qxy)
                scan.create_dataset('qz', data=self.qz)

                    # Add metadata
                scan.attrs['sample_name'] = self.sample_name
                scan.attrs['pi'] = self.Pi
                scan.attrs['scans'] = self.scans
                scan.attrs['energy'] = self.energy
                scan.attrs['alpha_i'] = self.alpha_i
                logger.info(f"2D image saved to {filename}")
        except Exception as e:
            logger.error(f"Scan already processed and saved to h5: {e}")

    def save_image_dat(self, filename=None):
        """
        Save the intensity map, qxy axis, and qz axis.
        """

        if filename is None:
            self._ensure_sample_dir()
            filename = self.saving_dir +'/GID_{}_scan_{}_2D.dat'.format(self.sample_name, self.scans)

        qxy_size = len(self.qxy)
        qz_size = len(self.qz)

        qxy_out = np.ravel(np.outer(self.qxy,np.ones(qz_size)))
        qz_out = np.ravel(np.tile(self.qz, qxy_size))
        image_out = np.ravel(self.data_gap)

        out = np.array([qxy_out, qz_out, image_out]).T
        try:
            np.savetxt(filename, out)
            logger.info(f"2D image saved to {filename}")
        except Exception as e:
            logger.error(f"Error saving GID image as .dat file: {e}")

    @staticmethod
    def line_step(x, slope, intercept, step_pos, step):
        line_y = intercept + x * slope
        step_y = np.zeros(len(x))
        step_y[np.argwhere(x > step_pos)] = 1 * step
        line_step = line_y + step_y
        return line_step

    @staticmethod
    def line_step_fitting(pars, x, data = None):
        # unpack parameters: extract .value attribute for each parameter
        parvals = pars.valuesdict()
        slope = parvals['slope']
        intercept = parvals['intercept']
        step = parvals['step']
        step_pos = parvals['step_pos']

        line_y = intercept + x * slope
        step_y = np.zeros(len(x))
        step_y[np.argwhere(x > step_pos)] = 1 * step
        model = line_y + step_y

        if data is None:
            return model
        return (model - data)

    @staticmethod
    def calibrate_mythen(filename, scanN, plot=True):

        detector_name = 'mythen2'
        angle_name = 'gam'

        logger.info('Calibrating Mythen from scan #{}'.format(scanN))
        try:
            with h5py.File(filename, "r") as f:
                # Using [()] to read dataset into numpy array immediately
                mythen = f.get(f"{scanN}.1/measurement/{detector_name}")[()]
                angle = f.get(f"{scanN}.1/measurement/{angle_name}")[()]

        except Exception as e:
            logger.error(f"Error loading scan {scanN}: {e}")
            raise

        beam_pos = np.argmax(mythen, axis=1)
        beam_int = np.max(mythen, axis=1)
        beam_pos = np.ma.masked_array(beam_pos, np.array(beam_int) < 2 * np.mean(mythen))

        gap_cen = (min(angle[~beam_pos.mask]) + max(angle[~beam_pos.mask])) / 2

        fit_params = lmfit.create_params(slope=-170, intercept=10000, step=100, step_pos=gap_cen)

        out = lmfit.minimize(GID.line_step_fitting, fit_params, args=(angle[~beam_pos.mask],),
                             kws={'data': beam_pos[~beam_pos.mask]})

        if plot:
            fig, ax = plt.subplots(1,1, figsize = (6,3), layout='tight')
            ax.plot(angle, beam_pos, 'o', alpha = 0.2)
            ax.plot(angle, GID.line_step_fitting(out.params, angle))

            ax.set_xlabel('gam, degree')
            ax.set_ylabel('beam position, px')
            ax.set_xlim(np.min(angle[~beam_pos.mask]), np.max(angle[~beam_pos.mask]))
            ax.set_ylim(np.min(beam_pos[~beam_pos.mask]), np.max(beam_pos[~beam_pos.mask]))

            plt.text(34, 500, f"PPD = {-out.params['slope'].value:.2f}, \nmythen_gap = {int(out.params['step'].value)}")
            plt.show()

        ppd = float(np.round(-out.params['slope'].value, 3))
        gap = int(out.params['step'].value)

        return ppd, gap

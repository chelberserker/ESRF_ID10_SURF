import copy
import logging
import os
import time
from math import sin, cos, pi
from typing import Optional, Tuple, Union, List, Dict, Any

import h5py
from datetime import datetime
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from ..base import BaseSurf

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

class GIXS(BaseSurf):
    """
    Class for processing GIS(W)AXS data obtained with a large 2D detector, e.g. Eiger, Pilatus
    """
    def __init__(self,
                 file,
                 scan,
                 alpha_i_name='chi',
                 detector_name='eiger4m',
                 monitor_name: str = 'mon',
                 transmission_name: str = 'autof_eh1_transm',
                 att_name: str = 'autof_eh1_curratt',
                 cnttime_name: str = 'sec',
                 energy_name='monoe',
                 I0 = 1e13,
                 saving_dir=None,
                 geometry = None,
                 **kwargs):
        """
        Initialize the GIXS processor.
        """
        super().__init__(file, scan, alpha_i_name=alpha_i_name, detector_name=detector_name,
                         monitor_name=monitor_name, transmission_name=transmission_name,
                         att_name=att_name, cnttime_name=cnttime_name, energy_name=energy_name,
                         I0=I0, saving_dir=saving_dir, **kwargs)

        self.geometry = geometry

        self.gi_integrator = None
        self.ip_range = None
        self.oop_range = None
        self.tilt_angle = None

        # Add GIXS specific appendable keys
        # Assuming GIXS scans might scan energy or angle, we append them.
        self.appendable_keys.extend(['data', 'alpha_i', 'energy'])

        self.load_data()

        # GIXS specific post-loading processing
        try:
            # Ensure alpha_i is an array of floats, possibly scalar wrapped in array if single value
            if np.size(self.alpha_i) == 1:
                self.alpha_i = np.array([float(self.alpha_i)])
        except Exception:
            logger.error('Error processing alpha_i')

    def _load_single_scan(self, scan_n:str):
        with h5py.File(self.file, "r") as f:
            base_path = f"{scan_n}.1"
            meas_path = f"{base_path}/measurement/"
            posi_path = f"{base_path}/instrument/positioners/"

            ### Trying to catch an exception here
            try:
                __incident_angle = f.get(f"{meas_path}{self.alpha_i_name}")
                if __incident_angle is None:
                    __incident_angle_path = posi_path
                else:
                    __incident_angle_path = meas_path
            except:
                # print(__incident_angle) # Logic from original code, likely debugging
                __incident_angle_path = posi_path

            data_dict = {
                'data': np.array(f.get(f"{meas_path}{self.detector_name}")),
                'alpha_i': np.array(f.get(f"{__incident_angle_path}{self.alpha_i_name}")),
                'monitor': np.array(f.get(f"{meas_path}{self.monitor_name}")),
                'transmission': np.array(f.get(f"{meas_path}{self.transmission_name}")),
                'attenuator': np.array(f.get(f"{meas_path}{self.att_name}")),
                'cnttime': np.array(f.get(f"{meas_path}{self.cnttime_name}")),
                'energy': np.array(f.get(f"{posi_path}{self.energy_name}")),
                'sample_name': str(f.get(f"{base_path}/sample/name/")[()])[2:-1:1],
            }


        logger.info('Loaded scan #%s', scan_n)
        return data_dict

    def __init_grazing__(self):
        ai = self.geometry
        self.gi_integrator = ai.promote(type_="pyFAI.integrator.fiber.FiberIntegrator")


    def integrate_single(self, frame=0, **kwargs):

        if self.gi_integrator is None:
            self.__init_grazing__()

        if len(self.alpha_i) == 1:
            __incident_angle = np.deg2rad(self.alpha_i)
        else:
            __incident_angle = np.deg2rad(self.alpha_i[frame])
        res_2d = self.gi_integrator.integrate2d_grazing_incidence(data=self.data[frame],
                                                                  sample_orientation=3,
                                                                  incident_angle=__incident_angle,
                                                                  normalization_factor = self.monitor[frame],
                                                                  **kwargs)

        self.ip_range = np.array([res_2d.inplane.min(), res_2d.inplane.max()])
        self.oop_range = np.array([res_2d.outofplane.max(),res_2d.outofplane.min()])

        return res_2d

    @staticmethod
    def __get_extent__(res_2d):
        ip_range = np.array([res_2d.inplane.min(), res_2d.inplane.max()])
        oop_range = np.array([res_2d.outofplane.max(), res_2d.outofplane.min()])
        extent = (ip_range[0], ip_range[1], oop_range[1], oop_range[0])
        return extent

    @staticmethod
    def plot_2d(res_2d, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(1,1,figsize=(10,10), layout='tight')

        extent = GIXS.__get_extent__(res_2d)
        ax.imshow(np.flipud(np.log10(res_2d[0])), extent=extent, **kwargs)

        ax.set_xlabel(res_2d.ip_unit.label)
        ax.set_ylabel(res_2d.oop_unit.label)


    def self_plot_2d(self, res_2d, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(1,1,figsize=(10,10), layout='tight')

        extent = self.__get_extent__(res_2d)
        ax.imshow(np.log10(res_2d[0]), extent=extent, **kwargs)



    def get_qz_cut(self, qxy_min, qxy_max, frame=0, **kwargs):
        """
        Extract a 1D cut along qxy by integrating over a range of qz.

        Args:
            qxy_min (float): Minimum qz value.
            qxy_max (float): Maximum qz value.
            frame (int): Frame number. Default is 0.

        Returns:
            integrate1d_grazing_incidence: Instance containing [qxy, intensity_profile].
                  intensity_profile is integrated intensity of counts in the qz range.
                  Contains all pyFAI metadata.
        """
        if self.gi_integrator is None:
            self.__init_grazing__()

        if len(self.alpha_i) == 1:
            __incident_angle = np.deg2rad(self.alpha_i)
        else:
            __incident_angle = np.deg2rad(self.alpha_i[frame])

        res_1d = self.gi_integrator.integrate1d_grazing_incidence(data=self.data[frame], incident_angle=__incident_angle,
                                                                 npt_oop=2000,
                                                                 npt_ip=200,
                                                                 sample_orientation=3,
                                                                 ip_range=(qxy_min, qxy_max),
                                                                 oop_range=self.oop_range,
                                                                 normalization_factor = self.monitor[frame],
                                                                  **kwargs)

        return res_1d

    def plot_qz_cut(self, qxy_min, qxy_max, frame=0, ax = None, plot_roi = False, **kwargs):
        res = self.get_qz_cut(qxy_min, qxy_max, frame, **kwargs)

        if plot_roi is False:
            if ax is None:
                fig, ax = plt.subplots(1,1,figsize=(6,6), layout='tight')

            ax.plot(res[0], res[1])
            ax.set_xlabel(res.unit.label)
            ax.set_ylabel('Intensity')

        else:
            if ax is not None:
                logger.error('plot_roi does not work with supplied axis!')
            else:
                fig, ax = plt.subplots(1,2,figsize=(12,6), layout='tight')


            res_2d_patch = self.integrate_single(frame, ip_range=(qxy_min, qxy_max))
            res_2d = self.integrate_single(frame)
            _v_min = np.log10(np.median(res_2d[0])-0.2*np.std(res_2d[0]))
            _v_max = np.log10(np.max(res_2d[0])+2*np.std(res_2d[0]))
            self.plot_2d(res_2d, ax=ax[0], vmin=_v_min, vmax=_v_max)
            self.plot_2d(res_2d_patch, ax=ax[0], vmin=_v_min, vmax=_v_max)

            ax[0].set_xlim(self.ip_range)
            ax[0].set_ylim(self.oop_range[1], self.oop_range[0])


            ax[1].plot(res[0], res[1])
            ax[1].set_xlabel(res.unit.label)
            ax[1].set_ylabel('Intensity')

        return fig, ax
